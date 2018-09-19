# fast-gmm

This package includes two components:
1. An implementation of ```CDDynamics```, as second order system which you can use to filter or generate smooth trajectories with the capability of setting a acceleration limits. 
How to use CDD-dynamics:  
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

How to use GMM:  
```
	Gaussians *GMM;
	GMM = new Gaussians(K, Dimention,"mu.txt","sigma.txt","prio.txt");  
	GMM->InitFastGaussians(0, 2); // I have no idea what is it!  
	start loop  
		Likelihood=GMM->GaussianProbFast(P)// P is a vector
	end loop  
```
How to use GMR:  
```
 	TODO
```
