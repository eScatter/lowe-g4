#!/bin/sh
job=results
for material in Alumina Silicon
do
	echo Starting job for ${material}
	sed "s/MATERIAL/${material}/g" YieldCurve.mac > Yield${material}.mac
	LD_PRELOAD="/usr/local/lib/libCADPhysics.so /usr/local/lib/libcsread.so /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5_cpp.so /home/annemarie/Downloads/hdf5-1.10.1/hdf5/lib/libhdf5.so" SEM YieldCurve.gdml Yield${material}.mac > Yield${material}.log 2>&1
	echo Done with job ${job} for ${material}
done
