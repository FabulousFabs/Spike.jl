using Spike
using Test

@testset "Spike.jl" begin
    @test Spike.greet_Spike() == "Hello!";
end
