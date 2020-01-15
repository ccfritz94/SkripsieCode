function[PLS_out]=PLS_output(Xcal,Ycal,Xtrain,Ytrain,Xtest,Ytest,nLVmax)
  
  dim_c=size(Xcal);
  dim_tr=size(Xtrain);
  dim_t=size(Xtest);
  
  PLS_out.B = plsN( Xcal,Ycal,nLVmax);
  
  %calculate model statistics
  Pred_train = (Xtrain)*PLS_out.B;
  
  
  Pred_test=(Xtest)*PLS_out.B;
 
RMSECV=[];
R2CV=[];
BIASCV=[];
SEPCV=[];
BIAS=[];
SEP=[];
RMSEP=[];
R2P=[];
RPDP=[];
RPDCV=[];
  
  for i = 1:nLVmax
      
    RMSECV(i)=(sum(((Pred_train(:,i)-Ytrain).^2)/(dim_tr(1)-1))^0.5);
    R2CV(i)=1-sum((Pred_train(:,i)-Ytrain).^2)/sum((Ytrain-mean(Ytrain)).^2);
    BIASCV(i)=sum(Pred_train(:,i)-Ytrain)/(dim_c(1));
    SEPCV(i)=(sum(((Pred_train(:,i)-Ytrain-BIASCV(i)).^2)/(dim_tr(1)-1)).^0.5);
    RPDCV(i)=std(Ytrain)/SEPCV(i);
    
    BIAS(i)=sum(Pred_test(:,i)-Ytest)/(dim_t(1));
    SEP(i)=(sum(((Pred_test(:,i)-Ytest-BIAS(i)).^2)/(dim_t(1)-1)).^0.5);
    RMSEP(i)=(sum(((Pred_test(:,i)-Ytest).^2)/(dim_t(1)-1)).^0.5);
    R2P(i)=1-sum((Pred_test(:,i)-Ytest).^2)/sum((Ytest-mean(Ytest)).^2);
    RPDP(i)=std(Ytest)/SEP(i);
  end
  
  [minCV,nLV]=min(RMSECV);
  
  PLS_out.stats=[nLV,RMSECV(nLV),R2CV(nLV),BIASCV(nLV),SEPCV(nLV),RPDCV(nLV),RMSEP(nLV),R2P(nLV),BIAS(nLV),SEP(nLV),RPDP(nLV)];
  PLS_out.RMSECV=RMSECV;
  PLS_out.Pred=Pred_test;
  PLS_out.RMSEP=RMSEP;
