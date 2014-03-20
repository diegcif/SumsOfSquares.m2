-- norm2.m2
-- 19.02.2004, Helfried Peyrl
-- Computes the Frobenius/Euclidean norm of a matrix
--
-- Input:  A:  matrix
-- Output: .:  squared Frobenius norm

norm2 = A -> (
     sqrt((sum apply(flatten entries A, i->i^2))_RR)
     )
