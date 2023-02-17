#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);
