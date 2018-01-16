from mpi4py import MPI
comm = MPI.COMM_WORLD
size=comm.Get_size()
rank = comm.Get_rank()
print("hello world from process {}".format(rank))
