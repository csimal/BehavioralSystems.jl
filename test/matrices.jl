
@testset "Matrices" begin
    @testset "Hankel Matrix" begin
        @test_throws DomainError hankel_matrix(1:5, 0)
        @test_throws DomainError hankel_matrix(1:5, -1)
        @test_throws DomainError hankel_matrix(1:3, 4)
        @test hankel_matrix(1:5, 1) == [1 2 3 4 5]
        @test hankel_matrix(1:5, 2) == [1 2 3 4; 2 3 4 5]
        @test hankel_matrix([1:5;], 2) == [1 2 3 4; 2 3 4 5]
        @test hankel_matrix([1:5;; 1:5]', 2) ==  [1 2 3 4; 1 2 3 4; 2 3 4 5; 2 3 4 5]
        @test hankel_matrix([1:5;; 1:5]', 3) == [1 2 3; 1 2 3; 2 3 4; 2 3 4; 3 4 5; 3 4 5]
        @test hankel_matrix([1:5,6:8], 2) == [hankel_matrix(1:5,2) hankel_matrix(6:8,2)]
    end
    @testset "antidiagonal" begin
        @test_throws DomainError antidiagonal(1,-1,1)
        @test_throws DomainError antidiagonal(1,1,-1)
        @test isempty(antidiagonal(0,2,4))
        @test isempty(antidiagonal(-1,2,4))
        @test isempty(antidiagonal(6,2,4))
        @test collect(antidiagonal(1,2,4)) == [(1,1)]
        @test collect(antidiagonal(2,2,4)) == [(1,2),(2,1)]
        @test collect(antidiagonal(3,2,4)) == [(1,3),(2,2)]
        @test collect(antidiagonal(5,2,4)) == [(2,4)]
    end
    @testset "Hankel Projection" begin
        local H_1 = hankel_matrix(1:5, 2)
        local H_2 = hankel_matrix([1:5,6:8], 2)
        @test_throws DomainError hankel_projection(H_1, 0)
        @test_throws DomainError hankel_projection(H_1, -1)
        @test hankel_projection(H_1,2) == float.(H_1)
        @test hankel_projection(H_1,1) == float.(H_1)
        @test hankel_projection(1:5, 1) == [1.0:5.0;;]
        @test hankel_projection([1 2 3 4 5], 1) == [1.0 2.0 3.0 4.0 5.0]
        @test hankel_projection(H_2, 2, [4,2]) == float.(H_2)
        @test hankel_projection(H_2, 1, [4,2]) == float.(H_2)
        @test hankel_projection(hankel_matrix([1:5;; 1:5]', 2), 2) == float.(hankel_matrix([1:5;; 1:5]', 2))
        @test hankel_projection(hankel_matrix([1:5;; 1:5]',3), 3) == float.(hankel_matrix([1:5;; 1:5]',3))
        @test hankel_projection(hankel_matrix([1:5;; 1:5]',3), 2) == [1.0 2.0 3.0; 1.0 2.5 3.5; 2.0 3.0 4.0; 2.0 3.0 4.0; 2.5 3.5 5.0; 3.0 4.0 5.0]
    end
    @testset "Multiplication Matrix" begin
        @test_throws DomainError multiplication_matrix(1:5, 4)
        @test_throws DomainError multiplication_matrix([1:5 6:10]', 4)
        @test multiplication_matrix(1:5, 5) == [1 2 3 4 5]
        @test multiplication_matrix(1:5, 6) == [1 2 3 4 5 0; 0 1 2 3 4 5]
        @test multiplication_matrix(1:5, 6, 1) == multiplication_matrix(1:5, 6)
        @test multiplication_matrix(1:4, 3, 2) == [1 2 3 4 0 0; 0 0 1 2 3 4]
        @test multiplication_matrix([1:5 6:10]', 6) == [1 2 3 4 5 0; 0 1 2 3 4 5; 6 7 8 9 10 0; 0 6 7 8 9 10]
        @test multiplication_matrix([1:4 5:8]', 3, 2) == [1 2 3 4 0 0; 0 0 1 2 3 4; 5 6 7 8 0 0; 0 0 5 6 7 8]
    end
end