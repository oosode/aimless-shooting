#!/bin/bash

Cqm=( $(grep Cqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))
Oqm=( $(grep Oqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))
Nqm=( $(grep Nqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))
Hqm=( $(grep Hqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))
Pqm=( $(grep Pqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))
Mgqm=( $(grep Mgqm qmmm.*.init | awk '{ for (i=3; i<=NF; i++) print $i};'))

mmkind_atoms=( $(grep MMKIND qmmm.*.init | awk '{print $2}'))
mmkind_value=( $(grep MMKIND qmmm.*.init | awk '{print $3}'))

link_qmid=( $(grep LINK qmmm.*.init | awk '{print $2}'))
link_kind=( $(grep LINK qmmm.*.init | awk '{print $3}'))
link_mmid=( $(grep LINK qmmm.*.init | awk '{print $4}'))
link_type=( $(grep LINK qmmm.*.init | awk '{print $5}'))
link_alph=( $(grep LINK qmmm.*.init | awk '{if (NF==6) {print $6} else {print "NULL"}}'))

e3x3_atoms=( $(grep G3X3 qmmm.*.init | awk '{print $2}'))
g3x3_exclu=( $(grep G3X3 qmmm.*.init | awk '{print $2}'))

echo ${mmkind_value[@]}
echo 
echo ${link_alph[@]}
echo ${Mgqm[@]}

echo ${#Hqm[@]}
