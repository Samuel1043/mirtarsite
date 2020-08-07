# run every param respectively

declare -a lr_arr=(0.001 0.01 0.1)
declare -a batch_arr=(64 32 16 4)
declare -a embed_arr=(100 200 50 400)
declare -a rnnlayer_arr=(1 2 3)
declare -a rnnhidden_arr=(100 200 50 400)
declare -a dropout_arr=(0 0.1 0.3 0.5)


for i in "${lr_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_lr$i/ --lr $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_lr$i'deepmirtar'/ --lr $i
done

for i in "${batch_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_batch$i/ --batch_size $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_batch$i'deepmirtar'/ --batch_size $i
done
for i in "${embed_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_embed$i/ --embedding_hidden $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_embed$i'deepmirtar'/ --embedding_hidden $i
done
for i in "${rnnlayer_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_rnnlayer$i/ --rnn_layer $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_rnnlayer$i'deepmirtar'/ --rnn_layer $i
done
for i in "${rnnhidden_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_rnnhidden$i/ --rnn_hidden $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_rnnhidden$i'deepmirtar'/  --rnn_hidden $i
done
for i in "${dropout_arr[@]}"; do
	python3 train_model.py ./dataset/ ./compare_dropout$i/ --class_dropout $i
	python3 train_model.py ./dataset_deepmirtar/ ./compare_dropout$i'deepmirtar'/ --class_dropout $i
done
