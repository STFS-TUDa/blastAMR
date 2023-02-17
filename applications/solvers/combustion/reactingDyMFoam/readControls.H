#include "readTimeControls.H"

correctPhi = pimple.dict().getOrDefault("correctPhi", true);

checkMeshCourantNo =
    pimple.dict().getOrDefault("checkMeshCourantNo", false);
