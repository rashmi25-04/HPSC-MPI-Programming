#Uncomment your choice of compiler
CC    =   mpicc
#CC    =   mpigcc
#CC    =   mpipgcc

FLAGS := 

all:	tree_based_broadcast	\
	tree_based_reduce	\
	allreduce_all_oneall	\
	butterfly_allreduce	\
	butterfly_allreduce_v1	\
	custom_allreduce_v3	\


clean:	rm tree_based_broadcast tree_based_reduce allreduce_all_oneall butterfly_allreduce butterfly_allreduce_v1 custom_allreduce_v

tree_based_broadcast: tree_based_broadcast.c
	$(CC) $(FLAGS) tree_based_broadcast.c -o tree_based_broadcast

tree_based_reduce: tree_based_reduce.c
	$(CC) $(FLAGS) tree_based_reduce.c -o tree_based_reduce

allreduce_all_oneall: allreduce_all_oneall.c
	$(CC) $(FLAGS) allreduce_all_oneall.c -o allreduce_all_oneall

butterfly_allreduce: butterfly_allreduce.c
	$(CC) $(FLAGS) butterfly_allreduce.c -o butterfly_allreduce

butterfly_allreduce_v1: butterfly_allreduce_v1.c
	$(CC) $(FLAGS) butterfly_allreduce_v1.c -o butterfly_allreduce_v1				

custom_allreduce_v3: custom_allreduce_v3.c
	$(CC) $(FLAGS) custom_allreduce_v3.c -o custom_allreduce_v3
