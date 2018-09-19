# fast-gmm

This package includes two components:
1. An implementation of ```CDDynamics```, a second order system which you can use to filter or generate smooth trajectories with the capability of setting a acceleration limits.  
Usage:  
```
	CDDynamics *testDyn;  
	testDyn = new CDDynamics(dim, dt, wn);  
	testDyn->SetVelocityLimits(velLimits);  
	testDyn->SetAccelLimits(AccelLimits);  
	testDyn->SetState(initial);  
	testDyn->SetTarget(target);  
	start loop  
		// if new target is set  
		testDyn->SetTarget(target);  
		// update dynamics  
		testDyn->Update();  
		// get state  
		testDyn->GetState(state);  
	end loop  
```

2. A fast implementation of sampling from GMM/GMR evaluation functions. Used for se-DS and lpv-DS motion generators.
Usage (GMM):
```
	Gaussians *GMM;
	GMM = new Gaussians(K, Dimention,"mu.txt","sigma.txt","prio.txt");  
	GMM->InitFastGaussians(0, 2); // I have no idea what is it!  
	start loop  
		Likelihood=GMM->GaussianProbFast(P)// P is a vector
	end loop  
```
Usage (GMR):
```
 	TODO
```

**NOTES:** This package also includes an implementation of SVR evaluation function and 3rd order polynomial trajectory generation. These might be useful for someone and in my opinion don't belong here. - Nadia
