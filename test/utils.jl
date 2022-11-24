
@testset "drss" begin
    
end
@testset "random_trajectory" begin
    local sys = drss(3,2,2)
    local T = 10
    local w = random_trajectory(sys, T)
    @test size(w,2) == T
    @test size(w,1) == 2+2
end
@testset "sizes" begin
    local sys = ss(rand(3,3), rand(3,2), rand(4,3), 0.0)
    @test sizes(sys) == (3,2,4)
end
@testset "canonical_permutation" begin
    @test canonical_permutation(1,1,3) == [1,4,2,5,3,6]
    @test canonical_permutation(0,1,3) == [1,2,3]
    @test canonical_permutation(1,0,3) == [1,2,3]
    @test isempty(canonical_permutation(0,0,3))
    @test isempty(canonical_permutation(1,1,0))
    @test canonical_permutation(2,2,2) == [1,2,5,6,3,4,7,8]
end