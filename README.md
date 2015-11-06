# fast-gmm

This package is an unofficial Seungsu's implementation of motion_generators package.   
There are several advantages of this package over the motion_generators package. First: It calculates GMM, GMR faster.   Second: CDD-dynamics is a second order system (You can set a limits of Acceleration)


How to use CDD-dynamics:  


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
  

How to use GMM:  


GMM = new Gaussians(K, Dimention,"mu.txt","sigma.txt","prio.txt");  
	GMM->InitFastGaussians(0, 2); // I have no idea what is it!  
start loop  
Likelihood=GMM->GaussianProbFast(P)// P is a vector3  
end loop  

