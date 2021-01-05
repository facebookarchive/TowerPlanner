%The limited version of MIDACO is licensed under a Creative Commons BY-NC-ND License [https://creativecommons.org/licenses/by-nc-nd/3.0/legalcode.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% GATEWAY HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           
%     _|      _|  _|_|_|  _|_|_|      _|_|      _|_|_|    _|_|    
%     _|_|  _|_|    _|    _|    _|  _|    _|  _|        _|    _|  
%     _|  _|  _|    _|    _|    _|  _|_|_|_|  _|        _|    _|  
%     _|      _|    _|    _|    _|  _|    _|  _|        _|    _|  
%     _|      _|  _|_|_|  _|_|_|    _|    _|    _|_|_|    _|_|  
%
%                                                   Version 6.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           See the MIDACO user manual for detailed information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Author (C) :   Dr. Martin Schlueter
%                   Information Initiative Center,
%                   Division of Large Scale Computing Systems,
%                   Hokkaido University, JAPAN.
%
%    Email :        info@midaco-solver.com
%
%    URL :          http://www.midaco-solver.com
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ solution ] = midaco( problem, option, key )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if 'option.parallel' field exists
if( isfield(option,'parallel') == 0 )
    option.parallel = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = option.parallel;
if(P <= 1)
    if( P == 0 )
        P = 1;
    end
    if( P < 0)
        P = -P;
        fprintf('\n\n*** SIMULATED Parallel Mode ***\n\n');   
    end
else
    %fprintf('\n*** Initializing MIDACO Parallel Mode ***\n');    
    if verLessThan('matlab','8.2')
      matlabpool close force local % Shut down existing parallel pool
      matlabpool % Start fresh parallel pool
    else 
      delete(gcp) % Shut down existing parallel pool   
      parpool % Start fresh parallel pool
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o    = problem.o;
n    = problem.n;
ni   = problem.ni;
m    = problem.m;
me   = problem.me;
xl   = problem.xl;
xu   = problem.xu;
x0   = problem.x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxeval     = option.maxeval;
maxtime     = option.maxtime;
printeval   = option.printeval;
save2file   = option.save2file;
param       = option.param;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-Scale Infinity-Values "Inf"
marker = 0;
for i = 1:n
  % Rescale XL
  if(isinf(xl(i)) == 1)
    xl(i)  = -1.0e+6;
    marker = 1; 
  end
  % Rescale XU 
  if(isinf(xu(i)) == 1)
    xu(i)  = 1.0e+6;
    marker = 1; 
  end 
  % Repair X0
  if(isinf(x0(i)) == 1)
    x0(i)  = xl(i);
    marker = 1; 
  end  
end
if(marker == 1)
  fprintf(1,'\n NOTE: Infinity Values have been rescaled to 1.0E+6 \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(1,n);
for i = 1:n
   x(i) = x0(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fff = zeros(1,P*o);
ggg = zeros(1,P*m);
xxx = zeros(1,P*n);
for c = 1:P
    for i = 1:n
        xxx((c-1)*n+i) = x0(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pfmax = 1000;
if( abs(param(10)) >= 1 ) 
  pfmax = abs(param(10));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iflag  = 0;
istop  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iout  = 1; % screen
if(save2file > 0)
    if(printeval > 0)
      iout1 = fopen('MIDACO_SCREEN.TXT','w+');
      iout2 = fopen('MIDACO_SOLUTION.TXT','w+');
    end
    if(save2file > 1)
    iout4 = fopen('MIDACO_HISTORY.TXT','w+');
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(printeval > 0)
mprint_head(P,o,n,ni,m,me,maxeval,maxtime,printeval,save2file,param,key,iout);
if(save2file > 0)
    
    mprint_head(P,o,n,ni,m,me,maxeval,maxtime,printeval,save2file,param,key,iout1);
    
    fprintf(iout2,' MIDACO - SOLUTION\n');
    fprintf(iout2,' -----------------\n');
    fprintf(iout2,'\n This file saves the current best solution X found by MIDACO.');
    fprintf(iout2,'\n This file is updated after every PRINTEVAL function evaluation,');
    fprintf(iout2,'\n if X has been improved.\n\n');
end
end
if(save2file > 1)
    mark = 1;
    hist_count = 0;
    save_history(mark,P,o,n,m,fff,ggg,xxx,iout4,istop);
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( P == 1 ) % Serial Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while(istop == 0) % Call MIDACO by Reverse Communication Loop
      
      % Evaluate Objectve Functions F(X) and Constraints G(X)
      [ f, g ] = problem.func( x );

      % Check and repain NaN and Inf values
      CHECK = isnan(f); for i=1:o; if (CHECK(i)==1); f(i)= 1.0e+99;end;end
      CHECK = isnan(g); for i=1:m; if (CHECK(i)==1); g(i)=-1.0e+99;end;end
      CHECK = isinf(f); for i=1:o; if (CHECK(i)==1); f(i)= 1.0e+88;end;end
      CHECK = isinf(g); for i=1:m; if (CHECK(i)==1); g(i)=-1.0e+88;end;end       
      
      % Call the MIDACO - MEX file
      [x,f,g,istop,iflag,evalx,timex,pc,px,pf,pg,pr,pp,acc,PARETOFRONT] =...
      midacox(P,o,n,ni,m,me,maxeval,maxtime,printeval,param,...
              xl,xu,key,istop,iflag,x,f,g);


      
      % MIDACO printing commands
      if(( pc > 0 )&&(printeval > 0))
          if(iflag >= 100); istop = 1; end        
          psize = PARETOFRONT(1);
          mprint_line(evalx, timex,pf,pr,iflag,o,psize,iout);

          if( o > 1 )
            for i=1:o
              dummyf(i) = pg(m+i);
              if( pg(m+o+i)>=1 )
                dummyf(i) = -pg(m+i);
              end
            end
            if(save2file>=1)
             printfront(PARETOFRONT,o,m,n,dummyf,pg,px,pfmax); 
            end
          end

          if(save2file > 0)
              mprint_line(evalx, timex,pf,pr,iflag,o,psize,iout1);
              if( (pc > 1)&&(iflag<0) )
                  fprintf(iout2,'\n\n            CURRENT BEST SOLUTION');
                  mprint_solution(o,n,m,me,evalx, timex, iflag,px,pf,pg,pr,pp,acc,psize,iout2);
              end            
              flush_output(iout1); 
              flush_output(iout2); 
          end
          flush_output(iout); 
      end
      if(save2file > 1)
        mark = 2;
        hist_count = hist_count + 1;
        if(hist_count>save2file) % ---> reset history
            mark = 1;
            hist_count = 0;
            fclose(iout4);
            iout4 = fopen('MIDACO_HISTORY.TXT','w+');
        end
        save_history(mark,P,o,n,m,f,g,x,iout4,istop);
      end         
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % Parallel Mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Create cells to be used within parfor-loop
  F_CELL = cell(P,o);
  G_CELL = cell(P,m);
  X_CELL = cell(P,n);

  blocks = 0; % counter for parallel processed blocks

  while(istop == 0) % Call MIDACO by Reverse Communication Loop
    
      blocks = blocks + 1;

      clear F_CELL;
      clear G_CELL;
      clear X_CELL;
            
      for c=1:P
          for i=1:n
              X_CELL{c}{i} = xxx( (c-1)*n+i );
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %           PARFOR - LOOP             %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      parfor c = 1:P                        %
                                            %
          x = zeros(1,n);                   %                         
          for i=1:n                         %
              x(i) = X_CELL{c}{i};          %
          end                               %
                                            %
          % Evaluate F(X) and  G(X)         %
          [ f , g ] = problem.func( x );    %  
                                            %
          for j=1:o                         %
              F_CELL{c}{j} = f(j);          %
          end                               %          
                                            %
          for j=1:m                         %
              G_CELL{c}{j} = g(j);          %
          end                               %
      end                                   %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for c=1:P
          for j=1:o
              fff( (c-1)*o+j ) = F_CELL{c}{j};
          end          
          for j=1:m
              ggg( (c-1)*m+j ) = G_CELL{c}{j};
          end
      end 

      % Check and repain NaN and Inf values
      CHECK = isnan(fff); for i=1:o*P; if (CHECK(i)==1); fff(i)= 1.0e+99;end;end
      CHECK = isnan(ggg); for i=1:m*P; if (CHECK(i)==1); ggg(i)=-1.0e+99;end;end
      CHECK = isinf(fff); for i=1:o*P; if (CHECK(i)==1); fff(i)= 1.0e+88;end;end
      CHECK = isinf(ggg); for i=1:m*P; if (CHECK(i)==1); ggg(i)=-1.0e+88;end;end       
       
      % Call the MIDACO - MEX file
      [xxx,fff,ggg,istop,iflag,evalx,timex,pc,px,pf,pg,pr,pp,acc,PARETOFRONT] =...
      midacox(P,o,n,ni,m,me,maxeval,maxtime,printeval,param,...
              xl,xu,key,istop,iflag,xxx,fff,ggg);
      
      % MIDACO printing commands
      if(( pc > 0 )&&(printeval > 0))
          if(iflag >= 100); istop = 1; end        
          psize = PARETOFRONT(1);
          mprint_line(evalx, timex,pf,pr,iflag,o,psize,iout);

          if( o > 1 )
            for i=1:o
              dummyf(i) = pg(m+i);
              if( pg(m+o+i)>=1 )
                dummyf(i) = -pg(m+i);
              end
            end
            if(save2file>=1)
             printfront(PARETOFRONT,o,m,n,dummyf,pg,px,pfmax); 
            end
          end

          if(save2file > 0)
              mprint_line(evalx, timex,pf,pr,iflag,o,psize,iout1);
              if( (pc > 1)&&(iflag<0) )
                  fprintf(iout2,'\n\n            CURRENT BEST SOLUTION');
                  mprint_solution(o,n,m,me,evalx, timex, iflag,px,pf,pg,pr,pp,acc,psize,iout2);
              end
              flush_output(iout1); 
              flush_output(iout2); 
          end
          flush_output(iout); 
      end
      if(save2file > 1)
        mark = 2;
        hist_count = hist_count + 1;
        if(hist_count>save2file) % ---> reset history
            mark = 1;
            hist_count = 0;
            fclose(iout4);
            iout4 = fopen('MIDACO_HISTORY.TXT','w+');
        end       
        save_history(mark,P,o,n,m,fff,ggg,xxx,iout4,istop);
      end       
  end
  
  x = xxx(1,1:n);
  g = ggg(1,1:m);
  f = fff(1,1:o);  

  % fprintf(iout,'\n ( Number of parallel processed  blocks: %i )\n',blocks);
  if(save2file > 0)
  % fprintf(iout1,'\n ( Number of parallel processed blocks: %i )\n',blocks);
  % fprintf(iout2,'\n ( Number of parallel processed blocks: %i )\n',blocks);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent test of solution given by MIDACO
% [ f, g ] = problem.func( x );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( o >= 2 )
    for i=1:2*o
        g(m+i) = pg(m+i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psize = PARETOFRONT(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(printeval > 0)
  if(o<=1)
      print_f = f;
  else
      print_f = pf;
  end  
  % Printing Final Information
  mprint_final(iflag,timex,evalx,maxtime,maxeval,param,iout)
  mprint_solution(o,n,m,me,evalx, timex, iflag,x,print_f,g,pr,pp,acc,psize,iout);
  if( save2file > 0)
      mprint_final(iflag,timex,evalx,maxtime,maxeval,param,iout1)
      mprint_solution(o,n,m,me,evalx, timex, iflag,x,print_f,g,pr,pp,acc,psize,iout1);
      mprint_final(iflag,timex,evalx,maxtime,maxeval,param,iout2)
      mprint_solution(o,n,m,me,evalx, timex, iflag,x,print_f,g,pr,pp,acc,psize,iout2);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      MIDACO Gateway Return Arguments to MAIN program
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solution.f = f;
solution.g = 0; %Dummy value in case no constraints are present
if(m>0)
  solution.g = g(1:m);
end
solution.x = x;
solution.time = timex;
solution.eval = evalx;
solution.iflag = iflag;
solution.psize = PARETOFRONT(1); %Number of pareto points in PARETOFRONT
solution.paretofront = PARETOFRONT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( printeval > 0 )
   fprintf(iout,'\n');  
   if( save2file > 0 )
       fclose all; % Disable this command, if problematic
       fprintf(iout,'\n Closing MIDACO Output Files (fclose all) \n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end % end of midaco main gateway code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   _____  _____  _____ _   _ _______ _____ _   _  _____   %%%%%%%%
%%%%%%  |  __ \|  __ \|_   _| \ | |__   __|_   _| \ | |/ ____|  %%%%%%%%
%%%%%%  | |__) | |__) | | | |  \| |  | |    | | |  \| | |  __   %%%%%%%%
%%%%%%  |  ___/|  _  /  | | | . ` |  | |    | | | . ` | | |_ |  %%%%%%%%
%%%%%%  | |    | | \ \ _| |_| |\  |  | |   _| |_| |\  | |__| |  %%%%%%%%
%%%%%%  |_|    |_|  \_\_____|_| \_|  |_|  |_____|_| \_|\_____|  %%%%%%%%
%%%%%%                                                          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mprint_head(P,o,n,ni,m,me,maxeval,maxtime,printeval,save2file,param,key,iout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* print head */
fprintf(iout,'\n MIDACO 6.0    (www.midaco-solver.com)');
fprintf(iout,'\n -------------------------------------\n\n');
fprintf(iout,' LICENSE-KEY:  %s',key);
fprintf(iout,'\n\n ----------------------------------------\n');
fprintf(iout,' | OBJECTIVES%5i | PARALLEL%10li |\n',o,P);
fprintf(iout,' |--------------------------------------|\n');
fprintf(iout,' | N%14li | MAXEVAL%11li |\n'  ,n,maxeval);
fprintf(iout,' | NI%13li | MAXTIME%11li |\n',ni,maxtime);
fprintf(iout,' | M%14li | PRINTEVAL%9li |\n',m,printeval);
fprintf(iout,' | ME%13li | SAVE2FILE%9li |\n',me,save2file);
fprintf(iout,' |--------------------------------------|\n');
dummy = 0;
for i=1:13
    if(param(i) ~= 0.0)
        dummy=1;
    end
end
if(dummy == 0)
    fprintf(iout,' | PARAMETER:    All by default (0)     |\n');
else
    if(param( 1)~=0.0); fprintf(iout,' | PARAM( 1) %14.7e ACCURACY    |\n',param( 1)); end
    if(param( 2)~=0.0); fprintf(iout,' | PARAM( 2) %14.7e RANDOM-SEED |\n',param( 2)); end
    if(param( 3)~=0.0); fprintf(iout,' | PARAM( 3) %14.7e FSTOP       |\n',param( 3)); end
    if(param( 4)~=0.0); fprintf(iout,' | PARAM( 4) %14.7e ALGOSTOP    |\n',param( 4)); end
    if(param( 5)~=0.0); fprintf(iout,' | PARAM( 5) %14.7e EVALSTOP    |\n',param( 5)); end
    if(param( 6)~=0.0); fprintf(iout,' | PARAM( 6) %14.7e FOCUS       |\n',param( 6)); end
    if(param( 7)~=0.0); fprintf(iout,' | PARAM( 7) %14.7e ANTS        |\n',param( 7)); end
    if(param( 8)~=0.0); fprintf(iout,' | PARAM( 8) %14.7e KERNEL      |\n',param( 8)); end
    if(param( 9)~=0.0); fprintf(iout,' | PARAM( 9) %14.7e ORACLE      |\n',param( 9)); end
    if(param(10)~=0.0); fprintf(iout,' | PARAM(10) %14.7e PARETOMAX   |\n',param(10)); end
    if(param(11)~=0.0); fprintf(iout,' | PARAM(11) %14.7e EPSILON     |\n',param(11)); end
    if(param(12)~=0.0); fprintf(iout,' | PARAM(12) %14.7e BALANCE     |\n',param(12)); end
    if(param(13)~=0.0); fprintf(iout,' | PARAM(12) %14.7e CHARACTER   |\n',param(13)); end        
end
fprintf(iout,' ----------------------------------------\n\n');

if( o == 1)
fprintf(iout,' [     EVAL,    TIME]        OBJECTIVE FUNCTION VALUE         VIOLATION OF G(X) \n');
fprintf(iout,' ------------------------------------------------------------------------------');
else
fprintf(iout,' [     EVAL,    TIME]   MULTI-OBJECTIVE PROGRESS   VIOLATION OF G(X)\n');
fprintf(iout,' -------------------------------------------------------------------   [PARETO]');
end
fprintf(iout,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mprint_line(evalx, timex,pf,pr,iflag,o,psize, iout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(iflag >= 10)
    if(iflag < 100)
        fprintf(iout,'\n *** WARNING ***   ( IFLAG =%4li ) \n\n',iflag);
    else
        fprintf(iout,'\n *** MIDACO INPUT ERROR ***   ( IFLAG =%4li ) \n',iflag);
        return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(abs(pf) <= 1.0e+10)
    if(abs(pr) <= 1.0e+5)
      if( o <= 1 ); fprintf(iout,' [%9i,%8i]        F(X):%19.8f         VIO:%13.6f\n',evalx, timex,pf,pr); end
      if( o >= 2 ); fprintf(iout,' [%9i,%8i]   PRO:%20.8f   VIO:%13.6f   [%6i]\n',evalx, timex,pf,pr,psize); end
    else
      if( o <= 1 ); fprintf(iout,' [%9i,%8i]        F(X):%19.8f         VIO:%13.6e\n',evalx, timex,pf,pr); end
      if( o >= 2 ); fprintf(iout,' [%9i,%8i]   PRO:%20.8f   VIO:%13.6e   [%6i]\n',evalx, timex,pf,pr,psize); end
    end
else
    if(abs(pr) <= 1.0e+5)
        if( o <= 1 ); fprintf(iout,' [%9i,%8i]        F(X):%19.8e         VIO:%13.6f\n',evalx, timex,pf,pr); end
      if( o >= 2 ); fprintf(iout,' [%9i,%8i]   PRO:%20.8e   VIO:%13.6f   [%6i]\n',evalx, timex,pf,pr,psize); end
    else
        if( o <= 1 ); fprintf(iout,' [%9i,%8i]        F(X):%19.8e         VIO:%13.6e\n',evalx, timex,pf,pr); end
      if( o >= 2 ); fprintf(iout,' [%9i,%8i]   PRO:%20.8e   VIO:%13.6e   [%6i]\n',evalx, timex,pf,pr,psize); end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mprint_solution(o,n,m,me,evalx, timex, iflag, px,pf,pg,pr,pp,acc,psize,iout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(iout,'\n --------------------------------------------\n');
fprintf(iout,' EVAL:%10li,  TIME:%8li,  IFLAG:%4li \n',evalx,timex,iflag);
fprintf(iout,' --------------------------------------------\n');

if( o <= 1 )
  if(abs(pf) <= 1.0e+18)
    fprintf(iout,' f(X) =%38.15f \n',pf);
  else
    fprintf(iout,' f(X) =%38.6e \n',pf);
  end
else

  if(abs(pf) <= 1.0e+18)
    fprintf(iout,' PROGRESS%36.15f \n',pf);
  else
    fprintf(iout,' PROGRESS%36.6e \n',pf);
  end

      fprintf(iout,' --------------------------------------------\n');
      fprintf(iout,' NUMBER OF PARETO POINTS%21i\n',psize);         
      fprintf(iout,' --------------------------------------------\n');   


  for i=1:o
    if(abs(pg(m+i)) <= 1.0e+18)
        if(pg(m+o+i)<1)
          fprintf(iout,' f(%4i) =%35.15f \n',i,pg(m+i));
        else
          fprintf(iout,' f(%4i) =%35.15f \n',i,-pg(m+i));
        end
    else
        if(pg(m+o+i)<1)
          fprintf(iout,' f(%4i) =%35.6e \n',i,pg(m+i));
        else
          fprintf(iout,' f(%4i) =%35.6e \n',i,-pg(m+i));
        end  
    end 

  end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(m > 0)
    
    fprintf(iout,' --------------------------------------------\n');
    
    if(iflag < 100)
        if(pr <= 1.0e+12)
            fprintf(iout,' VIOLATION OF G(X)%27.12f\n',pr);
        else
            fprintf(iout,' VIOLATION OF G(X)%27.6e\n',pr);
        end
        fprintf(iout,' --------------------------------------------\n');
    end
    
    for i=1:m
        if(i <= me)
            if(abs(pg(i)) <= acc)
                fprintf(iout,' g(%4i) =%16.8f  (EQUALITY CONSTR)\n',i,pg(i));
            else
                if(abs(pg(i)) <= 1.0e+7)
                    fprintf(iout,' g(%4i) =%16.8f  (EQUALITY CONSTR)  <---  INFEASIBLE  ( G NOT = 0 )\n',i,pg(i));
                else
                    fprintf(iout,' g(%4i) =%16.3e  (EQUALITY CONSTR)  <---  INFEASIBLE  ( G NOT = 0 )\n',i,pg(i));
                end
            end
        end
        if(i > me)
            if(pg(i) > -acc)
                if(abs(pg(i)) <= 1.0e+7)
                    fprintf(iout,' g(%4i) =%16.8f  (IN-EQUAL CONSTR)\n',i,pg(i));
                else
                    fprintf(iout,' g(%4i) =%16.3e  (IN-EQUAL CONSTR)\n',i,pg(i));
                end
            else
                if(abs(pg(i)) <= 1.0e+7)
                    fprintf(iout,' g(%4i) =%16.8f  (IN-EQUAL CONSTR)  <---  INFEASIBLE  ( G < 0 )\n',i,pg(i));
                else
                    fprintf(iout,' g(%4i) =%16.3e  (IN-EQUAL CONSTR)  <---  INFEASIBLE  ( G < 0 )\n',i,pg(i));
                end
            end
        end
    end
end
fprintf(iout,' --------------------------------------------         BOUNDS-PROFIL    \n');
for i=1:n
    if(abs(px(i)) <= 1.0e+14)
        if(pp(i) == 0);   fprintf(iout,' x(%4i) =%34.15f;  %%  XL___________________ \n',i,px(i)); end
        if(pp(i) == 1);   fprintf(iout,' x(%4i) =%34.15f;  %%  x____________________ \n',i,px(i)); end
        if(pp(i) == 2);   fprintf(iout,' x(%4i) =%34.15f;  %%  _x___________________ \n',i,px(i)); end
        if(pp(i) == 3);   fprintf(iout,' x(%4i) =%34.15f;  %%  __x__________________ \n',i,px(i)); end
        if(pp(i) == 4);   fprintf(iout,' x(%4i) =%34.15f;  %%  ___x_________________ \n',i,px(i)); end
        if(pp(i) == 5);   fprintf(iout,' x(%4i) =%34.15f;  %%  ____x________________ \n',i,px(i)); end
        if(pp(i) == 6);   fprintf(iout,' x(%4i) =%34.15f;  %%  _____x_______________ \n',i,px(i)); end
        if(pp(i) == 7);   fprintf(iout,' x(%4i) =%34.15f;  %%  ______x______________ \n',i,px(i)); end
        if(pp(i) == 8);   fprintf(iout,' x(%4i) =%34.15f;  %%  _______x_____________ \n',i,px(i)); end
        if(pp(i) == 9);   fprintf(iout,' x(%4i) =%34.15f;  %%  ________x____________ \n',i,px(i)); end
        if(pp(i) ==10);   fprintf(iout,' x(%4i) =%34.15f;  %%  _________x___________ \n',i,px(i)); end
        if(pp(i) ==11);   fprintf(iout,' x(%4i) =%34.15f;  %%  __________x__________ \n',i,px(i)); end
        if(pp(i) ==12);   fprintf(iout,' x(%4i) =%34.15f;  %%  ___________x_________ \n',i,px(i)); end
        if(pp(i) ==13);   fprintf(iout,' x(%4i) =%34.15f;  %%  ____________x________ \n',i,px(i)); end
        if(pp(i) ==14);   fprintf(iout,' x(%4i) =%34.15f;  %%  _____________x_______ \n',i,px(i)); end
        if(pp(i) ==15);   fprintf(iout,' x(%4i) =%34.15f;  %%  ______________x______ \n',i,px(i)); end
        if(pp(i) ==16);   fprintf(iout,' x(%4i) =%34.15f;  %%  _______________x_____ \n',i,px(i)); end
        if(pp(i) ==17);   fprintf(iout,' x(%4i) =%34.15f;  %%  ________________x____ \n',i,px(i)); end
        if(pp(i) ==18);   fprintf(iout,' x(%4i) =%34.15f;  %%  _________________x___ \n',i,px(i)); end
        if(pp(i) ==19);   fprintf(iout,' x(%4i) =%34.15f;  %%  __________________x__ \n',i,px(i)); end
        if(pp(i) ==20);   fprintf(iout,' x(%4i) =%34.15f;  %%  ___________________x_ \n',i,px(i)); end
        if(pp(i) ==21);   fprintf(iout,' x(%4i) =%34.15f;  %%  ____________________x \n',i,px(i)); end
        if(pp(i) ==22);   fprintf(iout,' x(%4i) =%34.15f;  %%  ___________________XU \n',i,px(i)); end
        if(pp(i) ==90);   fprintf(iout,' X(%4i) =%34.15f;  %%  WARNING: XL = XU      \n',i,px(i)); end
        if(pp(i) ==91);   fprintf(iout,' X(%4i) =%34.15f; ***ERROR*** (X > XU)      \n',i,px(i)); end
        if(pp(i) ==92);   fprintf(iout,' x(%4i) =%34.15f; ***ERROR*** (X < XL)      \n',i,px(i)); end
        if(pp(i) ==93);   fprintf(iout,' x(%4i) =%34.15f; ***ERROR*** (XL > XU)     \n',i,px(i)); end
        if(pp(i) < 0 );   fprintf(iout,' PROFIL-ERROR'); end
    else
        if(pp(i) == 0);   fprintf(iout,' x(%4i) =%34.1e;  %%  XL___________________ \n',i,px(i)); end
        if(pp(i) == 1);   fprintf(iout,' x(%4i) =%34.1e;  %%  x____________________ \n',i,px(i)); end
        if(pp(i) == 2);   fprintf(iout,' X(%4i) =%34.1e;  %%  _x___________________ \n',i,px(i)); end
        if(pp(i) == 3);   fprintf(iout,' x(%4i) =%34.1e;  %%  __x__________________ \n',i,px(i)); end
        if(pp(i) == 4);   fprintf(iout,' x(%4i) =%34.1e;  %%  ___x_________________ \n',i,px(i)); end
        if(pp(i) == 5);   fprintf(iout,' x(%4i) =%34.1e;  %%  ____x________________ \n',i,px(i)); end
        if(pp(i) == 6);   fprintf(iout,' x(%4i) =%34.1e;  %%  _____x_______________ \n',i,px(i)); end
        if(pp(i) == 7);   fprintf(iout,' x(%4i) =%34.1e;  %%  ______x______________ \n',i,px(i)); end
        if(pp(i) == 8);   fprintf(iout,' X(%4i) =%34.1e;  %%  _______x_____________ \n',i,px(i)); end
        if(pp(i) == 9);   fprintf(iout,' X(%4i) =%34.1e;  %%  ________x____________ \n',i,px(i)); end
        if(pp(i) ==10);   fprintf(iout,' X(%4i) =%34.1e;  %%  _________x___________ \n',i,px(i)); end
        if(pp(i) ==11);   fprintf(iout,' X(%4i) =%34.1e;  %%  __________x__________ \n',i,px(i)); end
        if(pp(i) ==12);   fprintf(iout,' X(%4i) =%34.1e;  %%  ___________x_________ \n',i,px(i)); end
        if(pp(i) ==13);   fprintf(iout,' X(%4i) =%34.1e;  %%  ____________x________ \n',i,px(i)); end
        if(pp(i) ==14);   fprintf(iout,' X(%4i) =%34.1e;  %%  _____________x_______ \n',i,px(i)); end
        if(pp(i) ==15);   fprintf(iout,' X(%4i) =%34.1e;  %%  ______________x______ \n',i,px(i)); end
        if(pp(i) ==16);   fprintf(iout,' x(%4i) =%34.1e;  %%  _______________x_____ \n',i,px(i)); end
        if(pp(i) ==17);   fprintf(iout,' X(%4i) =%34.1e;  %%  ________________x____ \n',i,px(i)); end
        if(pp(i) ==18);   fprintf(iout,' x(%4i) =%34.1e;  %%  _________________x___ \n',i,px(i)); end
        if(pp(i) ==19);   fprintf(iout,' X(%4i) =%34.1e;  %%  __________________x__ \n',i,px(i)); end
        if(pp(i) ==20);   fprintf(iout,' X(%4i) =%34.1e;  %%  ___________________x_ \n',i,px(i)); end
        if(pp(i) ==21);   fprintf(iout,' X(%4i) =%34.1e;  %%  ____________________x \n',i,px(i)); end
        if(pp(i) ==22);   fprintf(iout,' X(%4i) =%34.1e;  %%  ___________________XU \n',i,px(i)); end
        if(pp(i) ==90);   fprintf(iout,' X(%4i) =%34.1e;  %%  WARNING: XL = XU      \n',i,px(i)); end
        if(pp(i) ==91);   fprintf(iout,' X(%4i) =%34.1e; ***ERROR*** (X > XU)      \n',i,px(i)); end
        if(pp(i) ==92);   fprintf(iout,' X(%4i) =%34.1e; ***ERROR*** (X < XL)      \n',i,px(i)); end
        if(pp(i) ==93);   fprintf(iout,' X(%4i) =%34.1e; ***ERROR*** (XL > XU)     \n',i,px(i)); end
        if(pp(i) < 0 );   fprintf(iout,' PROFIL-ERROR'); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mprint_final(iflag,timex,evalx,maxtime,maxeval,param,iout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((iflag == 1)||(iflag == 2))
    
    if(timex >= maxtime); fprintf(iout,'\n OPTIMIZATION FINISHED  --->  MAXTIME REACHED'); end
    if(evalx >= maxeval); fprintf(iout,'\n OPTIMIZATION FINISHED  --->  MAXEVAL REACHED'); end
end
if((iflag == 3)||(iflag == 4))
    
    fprintf(iout,'\n OPTIMIZATION FINISHED  --->  ALGOSTOP (=%li)',param(4));
end
if((iflag == 5)||(iflag == 6))
    
    fprintf(iout,'\n OPTIMIZATION FINISHED  --->  EVALSTOP');
end
if((iflag == 7))
    
    fprintf(iout,'\n OPTIMIZATION FINISHED  --->  FSTOP REACHED');
end
fprintf(iout,'\n\n\n         BEST SOLUTION FOUND BY MIDACO');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printfront( PARETOFRONT,o,m,n,f,g,x,pfmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
ioutx = fopen('MIDACO_PARETOFRONT.tmp','w+');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(ioutx,'#########################################################\n');
fprintf(ioutx,'### This file contains the pareto front approximation ###\n');
fprintf(ioutx,'#########################################################\n');
fprintf(ioutx,'### Solution format:     F(1:O)    G(1:M)    X(1:N)   ###\n');
fprintf(ioutx,'#########################################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psize = PARETOFRONT(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(ioutx,'#                                       \n');
fprintf(ioutx,'#        O         M         N     PSIZE\n');
fprintf(ioutx,'#                                       \n');
fprintf(ioutx,'%10i%10i%10i%10i\n',o,m,n,psize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(ioutx,'#                       \n');
fprintf(ioutx,'#        MIDACO solution\n'); %%% DUMMY FOR NOW
fprintf(ioutx,'#                       \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:o
    dummy = f(i); xprint(dummy,ioutx);
end
for i=1:m
    dummy = g(i); xprint(dummy,ioutx);
end  
for i=1:n
    dummy = x(i); xprint(dummy,ioutx);
end   
fprintf(ioutx,'\n'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(ioutx,'#                       \n');
fprintf(ioutx,'#        All non-dominated solutions found by MIDACO\n');
fprintf(ioutx,'#                       \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:psize
    for i=1:o
        dummy = PARETOFRONT(2+o*(k-1)+i-1); 
        xprint(dummy,ioutx);
    end
    for i=1:m
        dummy = PARETOFRONT(2 + o*pfmax + m*(k-1)+i-1); 
        xprint(dummy,ioutx);
    end  
    for i=1:n
        dummy = PARETOFRONT(2 + o*pfmax + m*pfmax + n*(k-1)+i-1); 
        xprint(dummy,ioutx);
    end   
    fprintf(ioutx,'\n'); 
end
fclose(ioutx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rename file from temporary ".tmp" to text file ".TXT"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movefile('MIDACO_PARETOFRONT.tmp','MIDACO_PARETOFRONT.TXT','f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dummy=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_history( c,P,o,n,m,f,g,x,ioutx,istop) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent tix;
persistent previous_x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( c == 1 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tix = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
fprintf(ioutx,'###################################################\n');
fprintf(ioutx,'### This file contains the history of solutions ###\n');
fprintf(ioutx,'###################################################\n');
fprintf(ioutx,'### SOLUTION FORMAT:  F(1:O)  G(1:M)  X(1:N)    ###\n');
fprintf(ioutx,'###################################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(ioutx,'#                              \n');
fprintf(ioutx,'#        O         M         N \n');
fprintf(ioutx,'#                              \n');
fprintf(ioutx,'%10i%10i%10i\n',o,m,n);
fprintf(ioutx,'#                       \n');
fprintf(ioutx,'#        SOLUTION HISTORY (in chronological order)\n');
fprintf(ioutx,'#                       \n');  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( tix == 1 )
    tix  = 2;
  for i=1:o
    dummy = f(i); xprint(dummy,ioutx);
  end
  for i=1:m
    dummy = g(i); xprint(dummy,ioutx);
  end  
  for i=1:n
    dummy = previous_x(i); xprint(dummy,ioutx);
  end   

else

 for k=1:P
  for i=1:o
    if( istop <= 0)
      dummy = f((k-1)*o+i); xprint(dummy,ioutx);
    else
      dummy = f(i); xprint(dummy,ioutx);      
    end
  end
  for i=1:m
   if( istop <= 0)
      dummy = g((k-1)*m+i); xprint(dummy,ioutx);
    else
      dummy = g(i); xprint(dummy,ioutx);      
    end    
  end  
  for i=1:n
   if( istop <= 0)    
    dummy = previous_x((k-1)*n+i); xprint(dummy,ioutx);
   else
    dummy = x(i); xprint(dummy,ioutx);
   end
  end  
  
 end

end
fprintf(ioutx,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:P 
  for i=1:n
    previous_x((k-1)*n+i) = x((k-1)*n+i);
  end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flush_output(ioutx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xprint(dummy,ioutx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( abs(dummy) < 100000000 )
    fprintf(ioutx,'%19.7f',dummy);
else
    fprintf(ioutx,'%19.5e',dummy);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flush_output(unit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fflush(unit); % this command works only in Octave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

