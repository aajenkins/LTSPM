sffxcut = [np.arange(0,(scutcrop[1]-scutcrop[0])*dsize/dres,dsize/dres),ffdata[0][yconstant,scutcrop[0]:scutcrop[1]],ffdata[1][yconstant,scutcrop[0]:scutcrop[1]]]
sffycut = [np.arange(0,(scutcrop[1]-scutcrop[0])*dsize/dres,dsize/dres),ffdata[0][scutcrop[0]:scutcrop[1],xconstant],ffdata[1][scutcrop[0]:scutcrop[1],xconstant]]
