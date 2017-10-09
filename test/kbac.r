#!/usr/bin/Rscript

print(system.time(library("KBAC")))
print(system.time(a <- read.table("kbac.kbac",sep=" ")[,1:11]))
print(system.time(p <- KbacTest(a,alpha=9,num.permutation=1000)))
print(p)
#print(system.time(for (i in 1:10) { p <- KbacTest(a,alpha=2,num.permutation=1000) } ))
#print(system.time(for (i in 1:10) { a <- read.table("kbac.kbac",sep=" ")[,1:11] ; p <- KbacTest(a,alpha=2,num.permutation=1000) } ))
