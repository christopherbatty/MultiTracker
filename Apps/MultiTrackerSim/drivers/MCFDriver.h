//
//  MCFDriver.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MCFDriver__
#define __MCFDriver__

#include "surftrack.h"

class MCFDriver
{
public:
    static bool step(LosTopos::SurfTrack * st, double dt, bool bbwall = false);
    
private:
    static double determineMaxDt(const LosTopos::SurfTrack * st, const std::vector<LosTopos::Vec3d> & v);
    
    static void evaluateV(const LosTopos::SurfTrack * st, std::vector<LosTopos::Vec3d> & v);

    static int onBBWall(const LosTopos::Vec3d& pos);
    
};

#endif
