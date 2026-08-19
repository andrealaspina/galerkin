function ONE=MP()
% Unit of the working arithmetic: double(1) or, with multiprecision, mp(1)

try
  Digits=mp.Digits;
catch
  Digits=16;
end
if Digits==16
  ONE=1;
else
  ONE=mp(1);
end

end
