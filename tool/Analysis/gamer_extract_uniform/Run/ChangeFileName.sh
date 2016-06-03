#!/bin/bash

for target in Slice*; do new=SliceZ_$(echo $target | sed 's/^.*\(.\{6\}\)$/\1/'); mv $target $new; done
