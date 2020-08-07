import torch.nn as nn
import torch
import numpy as np


class SimilarityMatrixMask(nn.Module):
    def __init__(self, input_size=5, hidden_size=100, embedding_hidden=100, rnn_num_layer=1, classify_dropout=0, mir_len=30, site_len=50):
        super(SimilarityMatrixMask, self).__init__()

        def calculate_conv_maxpool_size(h, w, pad, dil, kernel, stride):
            hout = np.floor((h+2*pad[0]-dil[0]*(kernel[0]-1)-1)/stride[0]+1)
            wout = np.floor((w+2*pad[1]-dil[1]*(kernel[1]-1)-1)/stride[1]+1)
            return hout, wout

        # embedding layer for site and mir
        self.mirembedding = nn.Embedding(
            input_size, embedding_hidden, padding_idx=0).cuda()
        self.utrembedding = nn.Embedding(
            input_size, embedding_hidden, padding_idx=0).cuda()

        self.mirna = nn.LSTM(embedding_hidden, hidden_size,
                             num_layers=rnn_num_layer, batch_first=True, bidirectional=True)
        self.utr = nn.LSTM(embedding_hidden, hidden_size,
                           num_layers=rnn_num_layer, batch_first=True, bidirectional=True)

        # output layer for site and utr hidden state
        self.mirout = nn.Sequential(
            nn.Linear(hidden_size*2*rnn_num_layer, hidden_size), nn.ReLU())
        self.utrout = nn.Sequential(
            nn.Linear(hidden_size*2*rnn_num_layer, hidden_size), nn.ReLU())

        channels=3
        

        vertical_pad = (0, 0)
        vertical_dilation = (1, 1)
        vertical_kernal = (3, 30)
        vertical_stride = (1, 1)


        vertical_out = calculate_conv_maxpool_size(
            site_len, mir_len, vertical_pad, vertical_dilation, vertical_kernal, vertical_stride)

        ver_pool = (2, 1)

        vertical_out = calculate_conv_maxpool_size(
            vertical_out[0], vertical_out[1], vertical_pad, vertical_dilation, ver_pool, ver_pool)

        self.vertical = nn.Sequential(nn.Conv2d(1, channels, (3, 30), stride=(1, 1)), nn.BatchNorm2d(3),
                                      nn.ReLU(), nn.MaxPool2d((2, 1)), nn.Flatten())

        horizon_pad = (0, 0)
        horizon_dilation = (1, 1)
        horizon_kernal = (50, 3)
        horizon_stride = (1, 1)
        horizon_out = calculate_conv_maxpool_size(
            site_len, mir_len, horizon_pad, horizon_dilation, horizon_kernal, horizon_stride)

        hor_pool = (1, 2)

        horizon_out = calculate_conv_maxpool_size(
            horizon_out[0], horizon_out[1], horizon_pad, horizon_dilation, hor_pool, hor_pool)

        self.horizontal = nn.Sequential(nn.Conv2d(1, channels, (50, 3), stride=(1, 1)), nn.BatchNorm2d(3),
                                        nn.ReLU(), nn.MaxPool2d((1, 2)), nn.Flatten())

        flatten_feat = int(vertical_out[0] * \
            vertical_out[1]+horizon_out[0]*horizon_out[1])*3
        self.out = nn.Sequential(
            nn.ReLU(),
            nn.Dropout(classify_dropout),
            nn.Linear(hidden_size*2+flatten_feat, 1)
        )

    def forward(self, x):

        embed_utr = self.mirembedding(x[0].long().cuda())
        embed_mir = self.utrembedding(x[1].long().cuda())

        mir, (mir_hid, _) = self.mirna(embed_mir)
        utr, (utr_hid, _) = self.utr(embed_utr)

        utr_mask = x[2].long().cuda().unsqueeze(2).expand(utr.shape)
        mir_mask = x[3].long().cuda().unsqueeze(2).expand(mir.shape)

        # add kernal 1
        interaction = torch.bmm(
            utr*utr_mask, (mir*mir_mask).transpose(1, 2)).unsqueeze(1)
        interaction_ver = self.vertical(interaction)
        interaction_hor = self.horizontal(interaction)

        mir_hid = mir_hid.view(mir_hid.shape[1], -1)
        mir_hid = self.mirout(mir_hid)
        utr_hid = utr_hid.view(utr_hid.shape[1], -1)
        utr_hid = self.utrout(utr_hid)

        stacked = torch.cat(
            (mir_hid, utr_hid, interaction_ver, interaction_hor), 1)
        out = self.out(stacked)
        return out.squeeze(1)

    def saliency_interaction(self, x):
        '''
        plot saliency map 
        '''

        embed_utr = self.mirembedding(x[0].long().cuda())
        embed_mir = self.utrembedding(x[1].long().cuda())

        mir, (mir_hid, _) = self.mirna(embed_mir)
        utr, (utr_hid, _) = self.utr(embed_utr)

        utr_mask = x[2].long().cuda().unsqueeze(2).expand(utr.shape)
        mir_mask = x[3].long().cuda().unsqueeze(2).expand(mir.shape)

        # add kernal 1
        interaction = torch.bmm(
            utr*utr_mask, (mir*mir_mask).transpose(1, 2)).unsqueeze(1)
        interaction.retain_grad()
        
        interaction_ver = self.vertical(interaction)
        interaction_hor = self.horizontal(interaction)

        mir_hid = mir_hid.view(mir_hid.shape[1], -1)
        mir_hid = self.mirout(mir_hid)
        utr_hid = utr_hid.view(utr_hid.shape[1], -1)
        utr_hid = self.utrout(utr_hid)

        stacked = torch.cat(
            (mir_hid, utr_hid, interaction_ver, interaction_hor), 1)
        out = self.out(stacked)

        return out.squeeze(1), interaction, stacked
