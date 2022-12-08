
@testset "Behavior Representations" begin
    @testset "Model based Image Representation" begin
        A = [0 1; 1 0]; B = [2 0; 3 0]; C = [1 0]; D = [0 1]
        sys = ss(A,B,C,D,1)
        ℬ = ss2BT_modelbased(sys, nothing, 2)
        @test ℬ == [zeros(Int,4,2) I; C D zeros(Int,1,2); C*A C*A*B D]
    end
    @testset "Data-driven Image Representation" begin
        
    end
    @testset "Model based Kernel Representation" begin
        
    end
    @testset "Data-driven Kernel Representation" begin
        
    end
end