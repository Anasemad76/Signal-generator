clc;
clear all;
fs = str2double(input('ENTER SAMPLING FREQUENCY:', 's'));
while(fs<=0)
    fs=str2double(input('INVALID INPUT,PLEASE ENTER POSITIVE SAMPLING FREQUENCY:', 's'));
 end
    start_time = str2double(input('ENTER START TIME:', 's'));
      end_time = str2double(input('ENTER END TIME:', 's'));
      while(end_time<=start_time)
          end_time = str2double(input('INVALID INPUT!,Enter valid end time: ','s'));
      end
    breaknum = str2double(input('ENTER NUMBER OF BREAKPOINTS:', 's'));
    breakpoints = zeros(1, breaknum); % 1d array size of breaknum
    t_begin=start_time; % will help in breakpoints
     sign=zeros(1,breaknum+1); % array of signals to be drawn
for(i=1:breaknum)
          breakposition = str2double(input(['ENTER POSITION OF BREAKPOINT ' num2str(i)'], 's'));
        while(breakposition<=start_time || breakposition>=end_time)
                       breakposition= str2double(input(['INVALID BREAKPOINT, PLEASE ENTER POSITION OF BREAKPOINT ' num2str(i)'], 's'));
        end
         breakpoints(i)=breakposition; % position of breakpoints   
end
   breakpoints=sort(breakpoints); % after sorting breakpoints inorder
      

for(i=1:breaknum+1)
   
    fprintf("Signal Generator Menu:\n1-DC signal\n2-Ramp signal\n3-Polynomial Signal\n4-Exponential Signal\n5-Sinusoidal Signal\n");
    signal=input("Enter signal name",'s');
  
    while 1
       if(strcmpi(signal,'DC')==1||strcmpi(signal,'Ramp')==1||strcmpi(signal,'Polynomial')==1||strcmpi(signal,'Exponential')==1||strcmpi(signal,'Sinusoidal')==1)
        break;
       else
           signal=input("INVALID SIGNAL!,please enter correct signal name",'s');  
  
       end
       end
   

  if(strcmpi(signal,'DC')==1)
         sign(i)=1;
  else if(strcmpi(signal,'Ramp')==1)
            sign(i)=2;
  else if(strcmpi(signal,'Polynomial')==1)
            sign(i)=3;
   else if(strcmpi(signal,'Exponential')==1)
             sign(i)=4;
   else if(strcmpi(signal,'Sinusoidal')==1)
           sign(i)=5;






  end
  end
  end
  end
  end
end

  SignalArray=[];
 flagBreak=0; % to know if there is more breakpoints or not
 t_total=linspace(start_time,end_time,(end_time - start_time)*fs);
      PolyEq=0;
  for(i=1:length(breakpoints)+1)
  

      if(length(breakpoints)==0) % if no breakingpoint
          flagBreak=1; % no more breakpoints
          t_special1=start_time;
           t_special2=end_time;


      else if(i==length(breakpoints)+1) % if last breakingpoint
              flagBreak=1;
                t_special1=t_begin; % last breakpoint 
               t_special2=end_time; 
       
    
      else % general case in between signals
          flagBreak=0;
          t_special1=t_begin; %last breakpoint 0
          t_special2=breakpoints(i); % next breakpoint 2
         
      end
      end

          
      if(sign(i)==1) % if DC
           DCAmp=str2double(input('Enter DC Amplitude','s'));
           temp=DCAmp*ones(1,(t_special2-t_special1)*fs); 
           SignalArray=[SignalArray temp];
           if(flagBreak==0)
           t_begin=breakpoints(i);  
           end
         
      else if(sign(i)==2)% if Ramp
           RampSlope=str2double(input('Enter Ramp slope','s'));
           Ramp_Yinter=str2double(input('Enter Y intercept','s'));
           t_ramp=linspace(t_special1,t_special2,(t_special2-t_special1)*fs);
           temp=RampSlope*t_ramp +Ramp_Yinter;
            SignalArray=[SignalArray temp];
       if(flagBreak==0)
           t_begin=breakpoints(i);  
       end

      else if(sign(i)==3)%if poly
            t_Poly=linspace(t_special1,t_special2,(t_special2-t_special1)*fs);
         PolyPower=str2double(input('Enter Polynomial Power','s'));
           for(j=PolyPower:-1:1)
              PolyCoeff=str2double(input(['Enter Polynomial coefficient of ' num2str(j) ' power : '],'s'));
              PolyEq=(PolyEq+PolyCoeff.*(t_Poly).^j);
           end
             PolyAmplitude = str2double(input('Enter polynomial amplitude: ','s'));
            PolyInter=str2double(input('Enter Polynomial Intercept','s'));
           PolyEq=PolyAmplitude.*(PolyEq+PolyInter);
           temp=PolyEq;
            SignalArray=[SignalArray temp];
        if(flagBreak==0)
           t_begin=breakpoints(i);
           end
         


      else if(sign(i)==4)%if Exponential
           ExpAmp=str2double(input('Enter Exponential Amplitude','s'));
           ExpExpo=str2double(input('Enter Exponential Exponent','s'));
           t_Exp=linspace(t_special1,t_special2,(t_special2-t_special1)*fs);
          temp = ExpAmp*exp(ExpExpo*t_Exp);
           SignalArray=[SignalArray temp];
            if(flagBreak==0)
           t_begin=breakpoints(i);  
           end
          
      else if(sign(i)==5)%if sinusoidal
    
          while 1
              sinorcos=str2double(input('ENTER 1 for SIN AND 2 FOR COSINE:', 's'));
             if(sinorcos==1 || sinorcos==2)
                 break;
             end
          end
          SinuAmp=str2double(input('Enter Sinusiodal Amplitude','s'));
          SinuFreq=str2double(input('Enter Sinusoidal Freq ','s'));
          SinuPhase=str2double(input('Enter Sinusoidal Phase','s'));
          t_Sinu=linspace(t_special1,t_special2,(t_special2-t_special1)*fs);
          if(sinorcos==1)
            temp=SinuAmp*sin(2*pi*SinuFreq*t_Sinu+SinuPhase);
             SignalArray=[SignalArray temp];
               if(flagBreak==0)
           t_begin=breakpoints(i);  
           end
          else if(sinorcos==2)
              temp=SinuAmp*cos(2*pi*SinuFreq*t_Sinu+SinuPhase);
             SignalArray=[SignalArray temp];
          if(flagBreak==0)
           t_begin=breakpoints(i);  
           end
          

          end
          end
      end
      end
      end
      end
  end
  end
     
     


 plot(t_total,SignalArray); % draw signals 

 % part 2
SignalOperation=SignalArray;
t_Operation=t_total;
% acts as a stack to display each signal before modifing it
SignalOperationOld=SignalOperation; 
t_OperationOld=t_Operation;
while 1
 operation=str2double(input('Choose any operation to be done on the signal,if not press 7\n 1-Amplitude scaling\n 2-Time reversal\n 3-Time shift\n 4-Expanding\n 5-Compressing\n 7-Exit','s'));

while(operation<0 || operation>7)
    operation=str2double(input('INVALID INPUT!,please enter a number between 1 -> 6','s'));
 end

       subplot(2,2,[1,2]);
           plot(t_total,SignalArray);
           title('Original Signal');

 if(operation==1)
     Scalevalue=str2double(input('Enter Scaling value','s'));
     SignalOperation=Scalevalue.*SignalOperation;
     subplot(2,2,3);
     plot(t_OperationOld,SignalOperationOld);
     subplot(2,2,4);
     plot(t_Operation,SignalOperation);
     title('Scaling');

 else if(operation==2) 

         t_Operation=t_Operation.*-1;
              subplot(2,2,3);
     plot(t_OperationOld,SignalOperationOld);
           subplot(2,2,4);
         plot(t_Operation,SignalOperation)
         title('Time reversal');
         

 else if(operation==3) 
          Shiftvalue=str2double(input('Enter shift value','s'));
          t_Operation=t_Operation+Shiftvalue;
           subplot(2,2,3);
     plot(t_OperationOld,SignalOperationOld);
       subplot(2,2,4);
          plot(t_Operation,SignalOperation);
          title('Shifting');

   else if(operation==4)
           Expansion=str2double(input('Enter Expansion value','s'));
              t_Operation=t_Operation*Expansion;
               subplot(2,2,3);
     plot(t_OperationOld,SignalOperationOld);
       subplot(2,2,4);
              plot(t_Operation,SignalOperation);
              title('Expansion');

    else if(operation==5)
             Compression=str2double(input('Enter Compression value','s'));
                t_Operation=t_Operation/Compression;
                 subplot(2,2,3);
     plot(t_OperationOld,SignalOperationOld);
       subplot(2,2,4);
              plot(t_Operation,SignalOperation);
              title('Compression');

 

    else
     close all;
        disp("THANK YOU FOR USING SIGNAL GENERATOR")
        break;
    
    end


   end
 end
 end
 end
SignalOperationOld=SignalOperation; % to store the last signal before next opperation
t_OperationOld=t_Operation;
end