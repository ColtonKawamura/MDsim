export
    makeBox

function makeBox(boxX, boxY, cutoff, particleList)
    box = Box([boxX,boxY], cutoff)
    cl = CellList([p.position for p in particleList], box)
    return box, cl
end