
## ReadMe: Version 1.00

# A Static Power Grid State Estimator implemented as an Equivalent Circuit Programming (ECP) problem in MATLAB:


### What is an Equivalent Circuit Programming (ECP) problem?

ECP is a recently introduced framework for network optimization problems that was also applicable to power system simulations and optimizations when formulated using the equivalent circuit formulation (see [10]). The whole concept is based on representing the optimality conditions of a problem in terms of equivalent circuits and further solving the respective governing equations as a circuit simulation problem. What this means is that instead of using the traditional (existing) optimization algorithms, we applied the circuit simulation heuristics developed based on the physical characteristics of the optimization problem and its dual (called the adjoint network in circuit theory). This was demonstrated to enable scalable and efficient power grid optimizations with robust convergence properties. For more info on other application of ECP see [5]-[6].


### A new formulation for formulating the Static State Estimation Problem

The recently introduced Equivalent Circuit Formulation of defining the power flow problem based on current, voltage and admittance state variables [10],[6], demonstrated to enable the incorporation of more complex models within the power flow framework. For instance, the transient models such as the induction motor [7], or the more complex models of power electronic devices [8] can be now included within power flow without loss of generality. Most importantly, it was demonstrated that the measurement devices such as PMUs and RTUs can be simultaneously included within the same framework, which further enabled the new perspective on to static State Estimation problem. For more info see [1]-[4]


### Open Source Version of the Static State Estimator in MATLAB

We provide this basic version of the Static State Estimator for power grid as a supplemental material to the paper published and presented at the PowerTech Milan 2019 conference. For more information see [1].


 - System Requirements:
MATLAB version 7.3 (R2006b) or later

 - Downloading the code:
https://github.com/markojereminov/ECP_based_SE

 - Description of the code:
The provided code is based solely for State Estimation problem of balanced networks and has incorporated the 	following power grid models:
    1. PMU device
    2. RTU device
    3. In-Phase Transformer
    4. Phase Shifting Transformer
    5. Fixed Shunt device
    6. Transmission Line (PI-model)
    It is important to note that the provided code can handle both PMU and RTU within the same problem. For more 	information on the models see [1]. 

 - Test Cases: The input data for the State Estimation cases represents an .mpc file augmented with the 	additional structures related to PMU and RTU measurements. Both RTU and PMU structures are defined in terms 	of confidence intervals around the obtained measurements, i.e. PQV for RTUs and phasor voltage and current 	measurements for PMUs. See [1] for more info. Most importantly, the information about the actual measurement can be 	obtained by finding the mean of the confidence intervals. Lastly, we want to  emphasize that the test cases provided 	with the code do not represent the real measurement data but are synthetically created by running the power flow of the existing cases provided within MATPOWER [12] toolbox with the introduced noise and measurement 	errors as discussed in [1]-[2]


Running Simulations Syntax:

	1. EstimatedState = runSE(case6515rte_SE, NonlinearSEswitch); % this will run the linear SE problem [3]-[4] for  		the input case:case6515rte_SE if NonlinearSEswitch = 0,and it will solve a nonlinear version discussed in [1]-[2] 		if NonlinearSEswitch = 1.

	2. EstimatedState = runSE(casedata,NonlinearSEswitch,PMU_Gerr,RTU_weight,damp_coeff,epsilon,tol,maxIter); % in 			addition to the input case and SE switch, you can further set the other solver parameters, such as weights on the 		measurements (PMU_Gerr for PMU; and RTU_weight: for RTU), circuit diode heuristic parameters (damp_coeff), 				complementary slackness violation (diode coefficient - epsilon), and additional parameters for Newton Raphson solver 	(tol,maxIter). The more detailed description of the ‘runSE’ function is provided within the function itself 


 - Citing the code: If used for further research and comparisons please cite the following paper:
1. M. Jereminov, A. Jovicic, M. Wagner, G. Hug and L. Pileggi, “Equivalent Circuit Programming for Estimating the State of a Power System,” in Proc. IEEE PowerTech Milan, Milan, Italy, June 2019. 

### Disclaimer: 
The open source version of the code contains the most basic circuit simulation heuristic (diode limiting) which was demonstrated to enable robust, stable and efficient convergence properties for the all examined test cases. For more complex circuit simulation algorithms and homotopy methods see [9]-[11].

### License: 
This code is distributed as open-source under the terms provided in License file.

### Contact:
For more info about the code, models or issues contact the author:

    Marko Jereminov
    m.jereminov92@gmail.com
    Department of Electrical and Computer Engineering
    Carnegie Mellon University 
    Pittsburgh, PA 
    United States

References: 
 - [1] M. Jereminov, A. Jovicic, M. Wagner, G. Hug and L. Pileggi, “Equivalent Circuit Programming for Estimating the State of a Power System,” in Proc. IEEE PowerTech Milan, Milan, Italy, June 2019. 
 - [2] A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug An equivalent circuit formulation for power system state estimation including PMUs, 2018 North American Power Symposium (NAPS).
 - [3]  A. Jovicic, M. Jereminov, L. Pileggi, and G. Hug, A Linear Formulation for Power System State Estimation including RTU and PMU Measurements, in Proc. of ISGT Europe 2019, Sept. 2019, Bucharest, Romania.
 - [4] M. R. Wagner, M Jereminov, A. Pandey, L. Pileggi, "A Probabilistic Approach to Power System State Estimation using a Linear Algorithm," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy
 - [5] M. Jereminov, D. Bromberg, A. Pandey, M. Wagner, L. Pileggi, “Adjoint power flow analysis for evaluating feasibility”, IEEE Trans. on Smart Grid. (submitted).
 - [6] M. Jereminov, A. Pandey, L. Pileggi, "Equivalent Circuit Formulation for Solving AC Optimal Power Flow", in IEEE Transactions on Power Systems, 2019.
 - [7] A. Pandey, M. Jereminov, X. Li, G. Hug, L. Pileggi, “Unified Power System Analyses and Models using Equivalent Circuit Formulation, ”IEEE PES ISGT, Minneapolis, USA, 2016.
 - [8] M. Jereminov, A. Pandey, D. M. Bromberg, X. Li, G. Hug, L. Pileggi, “Steady-State Analysis of Power System Harmonics Using Equivalent Split-Circuit Models”, ISGT Europe, Ljubljana, October 2016.
 - [9] A. Pandey, M. Jereminov, M. R. Wagner, G. Hug and L. Pileggi, "Robust Convergence of Power Flow using Tx Stepping Method with Equivalent Circuit Formulation" in 2018 (PSCC), Dublin 2018.
 - [10] M. Jereminov, A Terzakis, M. R. Wagner, Amritanshu Pandey, Larry Pileggi, "Robust and Efficient Power Flow Convergence with G-min Stepping Homotopy Method," in Proc. of IEEE/EEEIC Genoa 2019, June 2019, Genoa, Italy.
 - [11] A. Pandey, M. Jereminov, M. R. Wagner, D. M. Bromberg, G. Hug- Glanzmann and L. Pileggi, "Robust Power Flow and Three Phase Power Flow Analyses," in IEEE Transactions on Power Systems. Vol. 34, pp. 616-626, 2019.
 - [12] R. Zimmerman, et.al., “MATPOWER: Steady-state operations, planning and analysis tools for power systems research and education”, IEEE Trans. on Power Systems, vol. 26, no.1, Feb 2011.
