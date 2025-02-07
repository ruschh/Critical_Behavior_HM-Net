Flavio R. Rusch, Antonio C. Roque, Osame Kinouchi

##**Influence of Topology on the Critical Behavior of Hierarchical Modular Neuronal Networks**


##Codes repository

This repository has source codes in the Fortran 95 language used to simulate hierarchical modular network topologies and their dynamics.

**ER_HM_NET.f95**: That code generates an HM network with neuron modules randomly and sparsely connected.

**KN_HM_NET.f95**: That code generates an HM network with neuron modules randomly connected with the same number of neighbors.

**FC_HM_NET.f95**: That code generates an HM network with neuron modules that are fully connected.

**NETWORK_DYNAMIC_ER_KN.f95**: That code simulates the network dynamics of HM-ER and KN-HM, calculating the quantities of interest.

**NETWORK_DYNAMIC_FC.f95**: That code simulates the network dynamics of FC-ER, calculating the quantities of interest.


## Running the code

To compile the code you must use the following command:

**f95 source_code.f95 -o exec.out**

After that, to run the code use the command: 

**./exec.out**

 


