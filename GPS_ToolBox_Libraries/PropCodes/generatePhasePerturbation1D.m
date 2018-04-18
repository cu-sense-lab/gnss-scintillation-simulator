function   [phase]=generatePhasePerturbation1D(nstep,ky,SF,thetap,StructureParams2D)
%Generate phase realization at nstep
global rootSDF1
    %Phase screen model
    if nstep==1  %Generate rootSDF as global variable
        Cp=StructureParams2D.Cp;
        p1=StructureParams2D.p1;
        q0=StructureParams2D.q0;
        p2=StructureParams2D.p2;
        rootSDF1=root_phaseSDF1(Cp,p1,p2,q0,ky,SF*sec(thetap)^2);
        phase=turbsim1(rootSDF1*sqrt(diff(ky(1:2))/2/pi));
    else
        phase=turbsim1(rootSDF1*sqrt(diff(ky(1:2))/2/pi));
    end
    return

