function Delay(DeltaT)
if(DeltaT>0) %end condition
    t=timer('timerfcn','Delay(0)','StartDelay',DeltaT);
    start(t);
    wait(t);
end