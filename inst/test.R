
devtools::load_all()
geSiteTable(.newdata=data.table(Gravels=c(90,95,100),Sands=c(5,10,25),Hs=c(90,100,110),Water=0,POP=100))


devtools::load_all()
getSiteProperties(Hs=c(100,50,131),USCS=c("GW","GP","GM","ML","SM"),NR=25,levels=c(0.16,0.50,"mean",0.84))

devtools::load_all()
getCylinderRoots(mo=0.45,lo=0.44)
getCylinderRoots(mo=0.45,lo=0.44001,model="lm")
getCylinderRoots(mo=0.45,lo=0.44001,model="nlm")
getCylinderRoots(mo=0.45,lo=0.44001,model="dt")
getCylinderRoots(mo=0.45,lo=0.44001,model="rf")
getCylinderRoots(mo=0.99,lo=0.51)
