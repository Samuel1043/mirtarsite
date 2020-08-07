# grid search (not used)
declare -a lr_arr=(0.001 0.01 0.1)
declare -a batch_arr=(32 16 4)
declare -a embed_arr=(100 200 50)
declare -a rnnlayer_arr=(1 2 3)
declare -a rnnhidden_arr=(100 200 50)
declare -a dropout_arr=(0 0.1 0.3)

for i in "${lr_arr[@]}"; do
    for i1 in "${batch_arr[@]}"; do
        for i2 in "${embed_arr[@]}"; do
            for i3 in "${rnnlayer_arr[@]}"; do
                for i4 in "${rnnhidden_arr[@]}"; do
                    for i5 in "${dropout_arr[@]}"; do
                        echo ./all_lr$i'_batch_'$i1'_embed_'$i2'_rnnlayer_'$i3'_rnnhidden_'$i4'_dropout_'$i5/
                    done
                done
            done
        done
    done
done
