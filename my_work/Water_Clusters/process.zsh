#!/bin/zsh

grep -E '(-*\d+\.\d+\s*\n*){6}' hydrogen.lammpstrj > H.dat
grep -E '(-*\d+\.\d+\s*\n*){6}' oxygen.lammpstrj > O.dat

