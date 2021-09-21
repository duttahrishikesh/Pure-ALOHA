# ALOHA and DMLLI-MAC Simulator

ALOHA.c:
This Code simuates Pure ALOHA behavior in MAC layer. It is an event-driven simulator where the simulation kernel is based on the concept of event scheduling, in terms of packet generation, transmission and reception. 

To run a simulation proceed using the following steps:

1. Input Number of Network nodes in No_of_Nodes
2. Define the network topology using the matrix topology[][], where topology[i][j]=1, if nodes i and j are connected. Else, topology[i][j]=0. To note, that the node index starts from 0.
3. Input the network traffic load (G) in Cap_G, if load distribution is homogeneous across all network nodes
4. In case of heterogeneous load, input individual node load g_i separately into small_gi, i=1,2,...,N
5. Observe the outputs in terms of networkwide throughput and individual node throughput

DMLLI-MAC.c:

To implement the proposed DMLLI MAC protocol, the Reinforcement Learning update equations are embedded on top of the pure ALOHA simulator. The procedure of execution is same as that of ALOHA.c.

The specific learning parameters learning rates, epsilon and discount factor can be defined by variables (alpha, beta), epsilon and gamma respectively
