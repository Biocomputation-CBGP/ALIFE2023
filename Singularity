BootStrap: docker
From: julia:1.8.5-alpine3.17

%post
	sudo apt-get update && sudo apt get -y install git
	julia -e 'using Pkg; Pkg.add(url="https://github.com/Biocomputation-CBGP/ALIFE2023.jl")'