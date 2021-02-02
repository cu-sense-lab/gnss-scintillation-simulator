function   [phase]=generatePhasePerturbation2D(nstep,Ky,Kz,SF,thetap,StructureParams)
%Generate phase realization at nstep
global rootSDF2
    %Phase screen model
    if nstep==1  %Generate rootSDF as global variable
        Cp=StructureParams.Cp;
        gnu1=StructureParams.gnu1;
        q0=StructureParams.q0;
        gnu2=StructureParams.gnu2;
        a=StructureParams.a;
        b=StructureParams.b;
        A=StructureParams.A;
        B=StructureParams.B;
        C=StructureParams.C;
        dky=diff(Ky(1,1:2));dkz=diff(Kz(1:2,1));
        rootSDF2=root_phaseSDF2(Cp,gnu1,gnu2,q0,Ky,Kz,SF*a*b*sec(thetap)^2,A,B,C);
        phase=turbsim2(rootSDF2*sqrt(dky/2/pi*dkz/2/pi));
    else
        dky=diff(Ky(1,1:2)); dkz=diff(Kz(1:2,1));
        phase=turbsim2(rootSDF2*sqrt(dky/2/pi*dkz/2/pi));
    end
    return

