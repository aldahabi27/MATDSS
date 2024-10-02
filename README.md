# MATDSS Application

MATDSS Application is a MATLAB-based tool that seamlessly integrates with [OpenDSS©](https://sourceforge.net/projects/electricdss/), a specialized standalone application for distribution network simulations. OpenDSS© provides an engine that interfaces with MATLAB via a COM interface, enabling efficient communication between the two platforms. To install OpenDSS, you can visit the [official OpenDSS page](https://sourceforge.net/projects/electricdss/) on SourceForge.

To run this application, you will need [MATLAB](https://www.mathworks.com/products/matlab.html), a high-level language and interactive environment for numerical computation, visualization, and programming. You can find more information or get a license from the [official MATLAB website](https://www.mathworks.com/products/matlab.html).
<!--
## Motivation

The development of MATDSS was driven by two main challenges:

1. **Complex Control Framework**: The control framework developed for this research resulted in a highly complex system with numerous parameters to manage. MATDSS addresses this by providing a graphical user interface (GUI) that allows for live monitoring, parameter management, and the ability to enable different components (such as low-pass filters, PD control, and disturbances). MATDSS also allows for modifying control parameters and plotting outputs for control areas and DERs directly within the application. A standout feature of MATDSS is its ability to partition a feeder without altering the simulation files of OpenDSS, thereby managing control actions and time-series simulations from within the same application.

2. **Lack of Time-Series Simulations in OpenDSS**: While OpenDSS is recognized as a standard solver for distribution systems, it lacks built-in time-series simulations for controlled structures like those proposed in this research. Integrating OpenDSS with MATLAB was essential to conduct these simulations and verify the presented work. This integration allows the framework to be easily adapted for different circuits and scales to control multiple areas within a feeder.
-->
## Key Features

MATDSS Application was developed to ease the integration of MATLAB with OpenDSS, especially for implementing complex control architectures. The development process focused on the following key features:

- **Modular Functions**: Allows extensive modification and overriding of default characteristics, enhancing adaptability for advanced control schemes.
- **Feeder Structure and Control Areas**: Supports partitioning feeders into multiple control areas and saving customized structures for easy access and modification.
- **DER and LC Configurations**: Provides comprehensive control over DER and LC configurations, including advanced control features like PID and low-pass filtering.
- **Simulation Setup Data**: Pre-defines simulation variables and parameters for quick recall, significantly reducing simulation setup time.
- **Additional Features**: Provides precise control over time configuration, real-time monitoring, customizable plotting, and export functions for comprehensive simulation analysis.

## GUI Overview

MATDSS offers an intuitive GUI that provides users the flexibility to switch between different OpenDSS circuits, configure simulations, and monitor outputs in real time. Below is a summary of the GUI components:

- **Main Window**: Lists available OpenDSS files, shows simulation plots, and provides configuration and OpenDSS file editor tabs.
- **Configuration Tabs**: Allows adding, editing, and saving DER configurations, managing disturbances, and setting up control areas.
- **Simulation Output Visualization**: Offers real-time monitoring of DER response, control areas' powers, voltages, currents, and associated dual variables with customizable plotting functionality.

## How to Use

To use MATDSS, follow these steps:

1. Load your OpenDSS feeder files into the application.
2. Configure DER, LC, and other control parameters through the GUI.
3. Set up simulation parameters and run simulations.
4. Monitor and visualize the simulation outputs in real-time and export the results as needed.

## Citing MATDSS Application

If you use the MATDSS Application in your research, please cite the following references:

1. I. Farhat, E. Ekomwenrenren, J. W. Simpson-Porco, E. Farantatos, M. Patel, and A. Haddadi, "A Multi-Area Architecture for Real-Time Feedback-Based Optimization of Distribution Grids," arXiv [eess.SY], 2024. Available: [https://arxiv.org/abs/2401.09694](https://arxiv.org/abs/2401.09694).

2. I. Farhat, MATDSS Application Ver 0.96, 2024. Available: [https://github.com/aldahabi27/MATDSS/tree/main](https://github.com/aldahabi27/MATDSS/tree/main).

3. I. Farhat, "Multi-Area Architecture for Real-Time Feedback-Based Optimization of Distribution Grids," Ph.D. dissertation, Dept. of Electrical and Computer Engineering, Univ. of Waterloo, Waterloo, ON, Canada, 2024. Available: [https://hdl.handle.net/10012/20848](https://hdl.handle.net/10012/20848).

## Contact

For any questions or support, please contact the developer at [ilyas.farhat@outlook.com](mailto:ilyas.farhat@outlook.com).

---

**MATDSS Application**  
Copyright (c) 2024, Ilyas Farhat  
This application is developed by 'Ilyas Farhat' as part of his research work towards his Ph.D. thesis at the University of Waterloo.


<meta name="google-site-verification" content="Z25mYn7q42VwL4sBivMvQs33wUQLQUtClwvwJPCPJe4" />
