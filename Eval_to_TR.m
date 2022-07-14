function TR = Eval_to_TR (NR, NU, Ntri,Eval,SNRT,TRT)

TR = zeros(Ntri,NR*NU);
for k=1:Ntri
    TR_MU = zeros(1,NR*NU);
    for n = 1: length(TRT)
        TR_MU(Eval(k,:) >= SNRT(n)) = TRT(n);
    end
    TR(k,:)=TR_MU;
end
%p126より、伝搬チャネル行列を変化させる各試行事に固有値のSNRを計算しており、この値がtr_tableのSNRの値を超えたらその変調方式を使用できるとして、TRTを使用できる%