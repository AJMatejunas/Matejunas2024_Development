function energyFE = func_scaleFEEnergy(energyFE,eFact)
    energyFE.artStrain = energyFE.artStrain*eFact;
    energyFE.damage = energyFE.damage*eFact;
    energyFE.extWork = energyFE.extWork*eFact;
    energyFE.internal = energyFE.internal*eFact;
    energyFE.penConstraint = energyFE.penConstraint*eFact;
    energyFE.penContact = energyFE.penContact*eFact;
    energyFE.kinetic = energyFE.kinetic*eFact;
    energyFE.strain = energyFE.strain *eFact;
    energyFE.totalBalance = energyFE.totalBalance*eFact;
    energyFE.viscDisp = energyFE.viscDisp*eFact;
end

