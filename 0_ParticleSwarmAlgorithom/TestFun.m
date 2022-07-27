function y = TestFun(X)

y = sin(X(1, :).^2 + X(2, :).^2);

end