#!/bin/sh
job=results
for material in Alumina Silicon
do
	echo Starting job for ${material}
	sed "s/MATERIAL/${material}/g" YieldCurve.mac > Yield${material}.mac
	LD_PRELOAD=/usr/local/lib/libCADPhysics.so SEM YieldCurve.gdml Yield${material}.mac > Yield${material}.log 2>&1
	echo Done with job ${job} for ${material}
done
