
@testset "drss" begin
    
end
@testset "random_trajectory" begin
    
end
@testset "sizes" begin
    local sys = ss(rand(3,3), rand(3,2), rand(4,3), 0.0)
    @test sizes(sys) == (3,2,4)
end