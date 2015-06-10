#include <stdio.h>
#include <AmberTopology.hpp>
#include <AmberRestart.hpp>
#include <AmberTrajectory.hpp>
#include <ErrorSystem.hpp>

int main(void)
{
    CAmberTopology     top;
    CAmberTrajectory   traj;
    CAmberRestart      crd;

    top.Load("test.parm7");
    crd.AssignTopology(&top);
    crd.Create();

    traj.AssignTopology(&top);
    traj.AssignRestart(&crd);

    traj.PrintInfo("test.netcdf",AMBER_TRAJ_UNKNOWN,AMBER_TRAJ_CXYZB,NULL);

    ErrorSystem.PrintErrors();

    return(0);
}
