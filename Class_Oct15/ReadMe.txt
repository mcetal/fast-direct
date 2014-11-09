Matlab Files from Class on Oct. 15

These files study the rank structure of the linear system arising from the exterior Laplace problem. 

Files for building linear system: 

ExtLaplace_MultConn_GGM.m: -solves the exterior Laplace problem outside 4 ellipses
                           -based on Greengard, Greenbaum, McFadden paper
                           
GetLinearSystemv2.m and GetLinearSystemv3.m: -supply the matrix to be compressed for the exterior Laplace problem.

MakeEllipse.m: -function used in creating linear system. 
               -contains contour information for each ellipse in the problem. 
                

Files for studying the ranks of the system: 

Ranks_SameDomain_ChangingN.m: -displays relationship between ranks of off-diagonal blocks and the number of points
                               on each ellipse

Ranks_ChangingDomain: -relationship between rank and bringing ellipses closer together

id_decomp.m: -function for ID from Fast Direct Solvers Workshop

Visualizing_SkeletonPoints.m: -displays skeleton points associated with the compression of horizontal and vertical blocks
                               
ID_Accuracy.m: -verifies the error in the ID with formulas given in Cheng et al. 


                             