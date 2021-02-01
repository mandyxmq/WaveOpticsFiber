
# WaveOpticsFiber

## Overview
This repository includes the major implementation of the SIGGRAPH Asia 2020 paper
**A Wave Optics Based Fiber Scattering Model**, by Mengqi (Mandy) Xia, Bruce Walter, Eric Michielssen, David Bindel and Steve Marschner.

More information about this project can be found at the [project website](https://mandyxmq.github.io/research/wavefiber.html).

We use the Boundary Element Method (BEM) to conduct a full wave simulation of plane wave scattered by an infinitely long cylinder (extruded from an arbitrary cross-section).

## Dependencies
This implementation depends on eigen3 (included as a submodule), OpenMP and the [complex_bessel library](https://github.com/joeydumont/complex_bessel). Please install OpenMP and and the complex_bessel library first.


## Usage
1. To compile the simulation code in Linux, first clone the repository
    ```
	git clone --recursive https://github.com/mandyxmq/WaveOpticsFiber.git
	cd WaveOpticsFiber
    ```

	Then download the [hankel function table](https://drive.google.com/drive/folders/1HWHP0rJZJ9aJb5pcStUWuWKUy9BktRcV?usp=sharing) and put the hankeldouble folder inside WaveOpticsFiber.

	Finally make a build directory and compile the code:
	```
	mkdir build
	cd build
	cmake ..
	make
	```

2. Run the simulation and generate scattering distribution, pdf, cdf, cross-section ratio data.<br/>
a. Elliptical cross-section
./Sim [name] [isEllipse] [semi-major radius (um)] [semi-minor radius (um)] [number of elements] [number of outgoing azimuthal angles] [number of incident longitudinal angles] [number of wavelength samples] [real part of the cylinder ior] [imaginary part of the cylinder ior]
	> **Note the direction conventions:** When the incident direction is perpendicular to the cylinder axis, the longitudinal angle is Pi/2; when the azimuthal angle between the incident and outgoing direction is Pi, it is the forward scattering direction.

	Examples:
	```
	./Sim circle 1 1 1 300 360 100 50 1.55 0
	```

	This command will simulate plane waves incident on a lossless circular cylinder of radius 1um, ior 1.55. It discretizes the cylinder perimeter using 300 line segments and simulates 50 wavelength samples raning from 400nm to 700nm and 100 incident theta angles ranging from 0 to Pi/2. There are 360 outgoing azimuthal angles ranging from 0 to 2Pi. The scattering distribution, pdf, cdf, cross-section ratio data will be saved to build/circle. 
	> **Note:** Due to the symmetry of the problem, only one incident azimuthal angle needs to be simulated and here the angle is 0.

	```
	./Sim ellipse 1 1.6 1 300 360 100 50 1.55 0
	```

	This command will simulate plane waves incident on a lossless elliptical cylinder. The semi-major radius is 1.6um, the semi-minor radius is 1um and the ior 1.55. It discretizes the cylinder perimeter using 300 line segments and simulates 50 wavelength samples raning from 400nm to 700nm and 100 incident theta angles ranging from 0 to Pi/2. There are 360 outgoing azimuthal angles ranging from 0 to 2Pi. The scattering distribution, pdf, cdf, cross-section ratio data will be saved to build/ellipse. 
	> **Note:** Due to the symmetry of the problem, the incident azimuthal angle ranges from 0 to Pi/2 and has 1/4 of the outgoing azimuthal samples. Thus there are 90 incident azimuthal angles in this case.

	<br/>
	b. Non-elliptical cross-section
	./Sim [name] [isEllipse] [x coordinate file] [y coordinate file] [number of outgoing azimuthal angles] [number of incident longitudinal angles] [number of wavelength samples] [real part of cylinder ior] [the file containing the imaginary part of the cylinder ior]

	Example:
	```
	 ./Sim silk 0 ../crosssection/xvec_silk.txt ../crosssection/yvec_silk.txt 360 100 50 ../ior/silkK_double.binary
	 ```
 
	This command will simulate plane waves incident on a lossy silk cylinder. The discretization information is contained in the files ../crosssection/xvec_silk.txt and ../crosssectionyvec_silk.txt. The simulation contains 50 wavelength samples raning from 400nm to 700nm, 100 incident theta angles ranging from 0 to Pi/2 and 360 outgoing azimuthal angles ranging from 0 to 2Pi. The scattering distribution, pdf, cdf, cross-section ratio data will be saved to build/silk. The incident azimuthal angle also ranges from 0 to 2Pi and has the same number of samples as the outgoing azimuthal samples. The wavelength dependent imaginary part of the cylinder ior is in ../ior/silkK_double.binary.

## Contact
If you have any questions regarding this work, please contact Mandy Xia (mandyxia AT cs.cornell.edu).
