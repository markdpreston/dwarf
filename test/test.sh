echo Plink Data
time dwarf plinkdata.dw
echo Plink 1
time plink --file plink --assoc --out plink &> /dev/null
echo Dwarf 0
time dwarf plinkPed0.dw &> /dev/null
echo Dwarf 10
time dwarf plinkPed10.dw 2> /dev/null
echo Dwarf 20
time dwarf plinkPed20.dw 2> /dev/null

echo KBAC Data
time dwarf kbacdata.dw
echo Dwarf 1
time dwarf kbac.dw
echo R 1
time ./kbac.r

echo Big test of all functionality: Null data
dwarf null.dw

echo Big test of all functionality: Rare data
dwarf rare.dw
