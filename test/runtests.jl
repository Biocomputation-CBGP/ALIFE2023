using Test
using ALIFE2023

@testset "Random genome generation -- genome is of proper length" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        n = rand(2nₒ:2:2048)
        genome = random_genome(N, n, 0, nₒ)

        @test length(genome) == n
    end
end

@testset "Random genome generation -- genome contains proper vertices" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        n = rand(2nₒ:2:2048)
        genome = random_genome(N, n, 0, nₒ)

        @test all(genome .<= N)
    end
end

@testset "Random genome generation -- no sources -- sinks are included in graph" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        n = rand(2nₒ:2:2048)
        genome = random_genome(N, n, 0, nₒ)

        @test all(x -> x ∈ genome, N-nₒ+1:N)
    end
end

@testset "Random genome generation -- no sources -- sinks are not used as inputs" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        n = rand(2nₒ:2:2048)
        genome = random_genome(N, n, 0, nₒ)

        @test all(x -> x ∉ genome[1:2:end], N-nₒ+1:N)
    end
end

@testset "Random genome generation -- all sources are included" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        nᵢ = rand(1:nₒ)
        n = rand(2*max(nₒ, nᵢ):2:2048)
        
        genome = random_genome(N, n, nᵢ, nₒ)
        @test all(x -> x ∈ genome, 1:nᵢ)

    end
end

@testset "Random genome generation -- no sources are outputs" begin
    for _ in 1:1024
        N = rand(2:128)
        nₒ = rand(1:N÷2)
        nᵢ = rand(1:nₒ)
        n = rand(2*max(nₒ, nᵢ):2:2048)
        
        genome = random_genome(N, n, nᵢ, nₒ)
        @test all(x -> x ∉ genome[2:2:end], 1:nᵢ)
    end
end

