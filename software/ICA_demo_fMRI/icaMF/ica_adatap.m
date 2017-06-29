function [S,A,loglikelihood,Sigma,chi,exitflag]=ica_adatap(X,prior,par,draw);

% ICA_ADATAP  Mean field independent component analysis (ICA)
%    [S,A,LL,SIGMA,CHI,EXITFLAG]=ICA_ADATAP(X) performs linear 
%    instanteneous mixing ICA by solving the adaptive TAP or mean field 
%    or variational mean field equations [1-4].
%   
%    Input and output arguments: 
%      
%    X        : Mixed signal
%    A        : Estimated mixing matrix
%    S        : Estimated source mean values                                   
%    SIGMA    : Estimated noise covariance
%    LL       : Log likelihood, log P(X|A,SIGMA), divided by  
%               number of samples - mean field estimate
%    CHI      : Estimated source covariances
%    EXITFLAG : 0, no convergence; 1,  maximum number of iterations  
%               reached and 2, convergence criteria met.  
%
%    Outputs are sorted in descending order according to 'energy'
%    sum(A.^2,1))'.*sum(S.^2,2).
%
%    [S,A,LL,SIGMA,CHI]=ICA_ADATAP(X,PRIOR) allows for specification 
%    of priors on A, S, SIGMA and the method to be used. PRIOR.method 
%    overrides other choices of priors. PRIOR.A, PRIOR.S, Prior.Sigma 
%    and PRIOR.method are character strings that can take the following 
%    values
%    
%    PRIOR.A :            constant      constant mixing matrix 
%                         free          free likelihood optimization
%                         positive      constrained likelihood 
%                                       optimization (uses QUADPROG)
%                         MacKay        corresponding to MacKay's 
%                                       hyperparameter update
%
%    PRIOR.S :            bigauss       sum of two Gaussians, default
%                                       variance 1 centered at +/-1.  
%                         binary        binary +/-1 
%                         binary_01     binary 0/1             
%                         combi         combinations of other priors
%                                       E.g. M1 binary_01's and M2
%                                       exponential's are specified by
%                                       PRIOR.S1='binary_01';
%                                       PRIOR.M1=M1; 
%                                       PRIOR.S2='exponential'; 
%                                       PRIOR.M2=M2;
%                         exponential   exponential (positive)
%                         heavy_tail    heavy_tailed (non-analytic 
%  				        power law tail)
%                         Laplace       Laplace (double exponential)
%                         Gauss         Gaussian for factor analysis
%                                       PRIOR.Sigma='diagonal' and
%                                       probabilistic PCA 
%                                       PRIOR.Sigma='isotropic'
%
%    PRIOR.Sigma :        constant      constant noise covariance
%                         free          free likelihood optimization 
%                         isotropic     isotropic noise covariance
%                                       (scalar noise variance)
%                         diagonal      different noise variance for 
%                                       each sensor.
%
%    PRIOR.method :       constant      constant A and Sigma 
%                                       for test sets. A and Sigma
%                                       should initialized, see 
%                                       below. Sets prior.A='constant'
%                                       and prior.Sigma='constant'.
%                         fa            factor analysis (FA) sets
%                                       prior.A='free', 
%                                       prior.S='Gauss' and
%                                       prior.Sigma='diagonal'.
%                         neg_kurtosis  negative kurtosis set
%                                       prior.S='bigauss'.
%                         pos_kurtosis  positive kurtosis sets
%                                       prior.S='heavy_tail'. 
%                         positive      positive ICA sets
%                                       prior.S='exponential' and
%                                       prior.A='positive'.
%                         ppca          probabilistic PCA (pPCA)
%                                       sets prior.A='free', 
%                                       prior.S='Gauss' and 
%                                       prior.Sigma='isotropic'.
%
%    Default values are PRIOR.A='free', PRIOR.S='Laplace',
%    PRIOR.Sigma='isotropic' and PRIOR.method not set.
%
%    [S,A,LL,SIGMA,CHI]=ICA_ADATAP(X,PRIOR,PAR) is used for giving 
%    additional arguments. PAR is a structure array with the 
%    following fields 
%   
%       sources           Number of sources (default is quadratic
%                         size(X,1)). This parameter is overridden 
%                         by size(PAR.A_init,2) if PAR.A_init is
%                         defined.
%       solver            method for solving mean field equations
%                         can take the following character string
%                         values, default is sequential: 
%
%                         beliefprop2   Belief propagation -      
%                                       Sequential second order 
%                                       method by T. Minka. 
%                         sequential    normal sequential iterative
%                                       update (first order method) 
%                         
%       A_init            Initial mixing matrix (default is Toeplitz
%                         type with small randmom values).        
%       S_init            Initial source values (default is zero).
%       Sigma_init        Initial noise covariance (scalar for 
%                         isotropic). Default is Sigma_rel times
%                         empirical covariance.
%       Sigma_rel         See above, default 1 for PRIOR.A=
%                         'A_positive' and 100 otherwise. 
%       max_ite           Maximum of A and SIGMA updates (default
%                         is 20)
%       S_max_ite         Maximum number of S update steps in each
%                         E-step (default is 100)
%       tol               Termination tolerance for program in 
%                         terms of relative change in 1-norm of 
%                         A and Sigma (default is 10^-3).
%       S_tol             Termination tolerance for S updates in 
%                         terms of the squared error of the S
%                         mean field equations (default is 10^-10).
%
%    [S,A,LL,SIGMA,CHI]=ICA_ADATAP(X,PRIOR,PAR,DRAW) with DRAW 
%    different from zero will output runtime information. DRAW
%    equal to 2 will output the runtime value of the log 
%    likelihood. Default value is 1.
%
%    The memory requirements are O(M^2*N+D*N), where D is the 
%    number of sensors, N the number of samples (X is D*N) and M is 
%    the number of sources. The computational complexity is 
%    O(M^3*N) plus for positive mixing matrix, D times the 
%    computational complexity of a M-dimensional quadratic 
%    programming problem.   
%
%    It is recommended as a preprocessing step to normalize the 
%    data, e.g. by dividing by the largest element X=X/max(max(X)).

%
% References [1-4]:
%
%    Tractable Approximations for Probabilistic Models: The Adaptive 
%    Thouless-Anderson-Palmer Mean Field Approach
%    M. Opper and O. Winther
%    Phys. Rev. Lett. 86, 3695-3699 (2001).
%
%    Adaptive and Self-averaging Thouless-Anderson-Palmer Mean Field 
%    Theory for Probabilistic Modeling
%    M. Opper and O. Winther
%    Physical Review E 64, 056131 (2001). 
% 
%    Mean Field Approaches to Independent Component Analysis
%    Pedro A.d.F.R. Højen-Sørensen, Ole Winther and Lars Kai Hansen
%    Neural Computation 14, 889-918 (2002).
%  
%    TAP Gibbs Free Energy, Belief Propagation and Sparsity
%    L. Csato, M. Opper and O. Winther
%    In Advances in Neural Information Processing Systems 14 (NIPS'2001), 
%    MIT Press (2002).
 
% - by Ole Winther 2002 - IMM, Technical University of Denmark
% - http://isp.imm.dtu.dk/staff/winther/
% - version 2.2

% Uses adaptive TAP convention refs. [1,2] besides in the definition of the mean 
% function where it uses the ICA-convention of ref. [3]. 

[D,N]=size(X);                                                % number of sensor, samples
try prior=prior; catch prior=[]; end                          % make sure prior is defined
%try par=par; catch par=[]; end                               % make sure par is defined
[prior]=method_init(prior);                                   % initialize method
[A,A_cmd,A_prior,M]=A_init(prior,par,D);                      % initialize A
A_norm_old=norm(abs(A),1);                                    % quantify chance in A
[S,S_mean_cmd_all,S_mean_cmd_seq,likelihood_cmd,likelihood_arg,...
 S_arg_seq,cS_prior,cM,S_prior]=S_init(prior,par,M,N);        % initialize S
[Sigma,Sigma_cmd,Sigma_eta,Sigma_prior,Sigma_rel]=...
Sigma_init(prior,par,A_prior,S_prior,X,A,S);                  % initialize Sigma
dim_Sigma=size(Sigma,1)*size(Sigma,2);                        % number of elements in Sigma 
Sigma_norm_old=norm(abs(Sigma),1);                            % quantify chance in Sigma.
try max_ite=par.max_ite; catch max_ite=20; end                % default maximum number of EM-steps 
try S_max_ite=par.S_max_ite; catch S_max_ite=100; end         % default maximum number of S updates in each E-step 
try S_min_ite=par.S_min_ite; catch S_min_ite=1; end           % default minimum number of S updates in each E-step 
try tol=par.tol; catch tol=10^-3; end                         % termination tolerance relative chance in norm(Sigma)+norm(A) 
try S_tol=par.S_tol; catch S_tol=10^-10; end; S_tolN=S_tol/N; % termination tolerance chance in dS.*dS 
try                                                           % default solver is sequential
  switch par.solver
    case {'beliefprop2','convergent'}
      solver=par.solver;    
    otherwise
      solver='sequential';
  end
catch solver='sequential'; end
try draw=draw; catch draw=1; end                              % set draw 

% here starts the algorithm! 

if strcmp(solver,'beliefprop2')
  fprintf(' ICA Adaptive TAP - %s solver\n',solver); 
else  
  fprintf(' ICA Linear Response - %s solver\n',solver); 
end
fprintf(' %s source prior\n',S_prior); %s A and %s noise optimization\n',S_prior,A_prior,Sigma_prior); 
if strcmp(S_prior,'combi')
  for i=1:length(cM) 
    if cM(i)==1 fprintf('   1 %s source prior\n',cS_prior{i});  
    else fprintf('   %i %s source priors\n',cM(i),cS_prior{i}); end
  end
end
fprintf(' %s A and %s noise optimization\n',A_prior,Sigma_prior); 
fprintf(' sensors: D = %i\n sources: M = %i\n samples: N = %i\n\n',D,M,N); 

dS=zeros(M,N); S_old=zeros(M,N);
chi=zeros(M,M,N); diagchi=zeros(M,N);
G=zeros(M,M,N);
diagG=zeros(M,N);
hpred=zeros(M,1); V=zeros(M,N); dV=zeros(M,N); 
traceSS=zeros(M,M); tracechi=zeros(M,M);

dA_rel=Inf; dSigma_rel=Inf;

if draw 
  if dim_Sigma==D
    fprintf(' noise variance = %g \n\n',sum(Sigma)/D);
  else       
    fprintf(' noise variance = %g \n\n',trace(Sigma)/size(Sigma,1));
  end  
end

EM_step=0;

while EM_step<max_ite & (dA_rel>tol | dSigma_rel>tol) 

  EM_step=EM_step+1;

  if dim_Sigma==D
    InvSigma=1./Sigma;  
    for i=1:M  % set up matrices for M-step.
      for j=i:M J(i,j)=-(A(:,i))'*(InvSigma.*A(:,j)); J(j,i)=J(i,j); end
      for k=1:N h(i,k)=(A(:,i))'*(InvSigma.*X(:,k)); end
    end      
  else
    InvSigma=inv(Sigma); 
    J=-A'*InvSigma*A;
    h=A'*InvSigma*X;
  end

  if EM_step>1 % heuristic for initializing V so that V<V_0
    Vrel=repmat(diag(J)./diagJ,1,N); V=Vrel.*V;
  else 
    V=zeros(M,N);
  end
  diagJ=diag(J);  

  V_0=-repmat(diag(J),1,N);
  J=J-diag(diag(J)); 

  S_ite=0;  
  dSdS=Inf;
  dSdSN=Inf*ones(N,1);
  I=1:N;

  switch solver 

    case 'beliefprop2'

      alpha=S;     
      Omega=zeros(M,N);
      for i=1:N
        G(:,:,i)=J;
      end

      while (~isempty(I) & dSdS>S_tol & S_ite<S_max_ite)

        S_ite=S_ite+1;
  
        indx=randperm(M);

        for sindx=1:M

          cindx=indx(sindx); % current index 

          % update hcav and V
          if M==1
            hpred=(squeeze(G(1,1,I)))'.*alpha(1,I);
          else
            hpred=sum(shiftdim(G(cindx,:,I)).*alpha(:,I),1);
          end
          Gc=squeeze(G(cindx,cindx,I))';
          dV(cindx,I)=Gc./(1+Gc.*Omega(cindx,I))-V(cindx,I);
          if S_ite>1
            Vfrac=0.7; Vind=(abs(dV(cindx,I))<=Vfrac*(V_0(cindx,I)-V(cindx,I)));
            dV(cindx,I)=Vind.*dV(cindx,I)+Vfrac*(1-Vind).*sign(dV(cindx,I)).*(V_0(cindx,I)-V(cindx,I)); 
            V(cindx,I)=V(cindx,I)+dV(cindx,I);
          end 
          hcav=hpred-V(cindx,I).*(Omega(cindx,I).*hpred+alpha(cindx,I)); 
       
          % update S and derivative
          S_old=S(cindx,I);
          eval(S_mean_cmd_seq); S(cindx,I)=f;  % e.g [f,dfdg] = exponential(h(cindx,I)+hcav,V_0(cindx,I)-V(cindx,I));
          dS(cindx,I)=S(cindx,I)-S_old;

          % update alpha and Omega
          alpha_old=alpha(cindx,I); alpha(cindx,I)=(S(cindx,I)-dfdg.*hcav)./(1+dfdg.*V(cindx,I));
          Omega_old=Omega(cindx,I); Omega(cindx,I)=dfdg./(1+dfdg.*V(cindx,I));                   
          
          % update G 
          for j=1:length(I)
            i=I(j); dOmega=Omega(cindx,i)-Omega_old(j); sG=squeeze(G(:,:,i));
            if any((diag(sG)+dOmega/(1-dOmega*sG(cindx,cindx)-eps)*sG(:,cindx).^2)<zeros(M,1)) 
              dOmega=0; Omega(cindx,i)=Omega_old(j)+dOmega;
              alpha(cindx,i)= 0.5*(alpha_old(j)+alpha(cindx,i));
            end
            if (dOmega~=0) % use Sherman-Woodbury
              b = sG(:,cindx); G(:,:,i)=sG+dOmega/(1-dOmega*sG(cindx,cindx))*b*b'; 
            end
          end
 
        end % over variables - sindx

        dSdSN=sum(dS.*dS,1);
        I=find(dSdSN > S_tolN);
        dSdS=sum(dSdSN);

      end; % epoch loop - S_ite 
 
    case 'sequential'
     
      while (~isempty(I) & dSdS>S_tol & S_ite<S_max_ite)     

        S_ite=S_ite+1;  

        indx=randperm(M);

        for sindx=1:M

          cindx=indx(sindx); % current index 

          % update hcav           
          hcav=J(cindx,:)*S(:,I)-V(cindx,I).*S(cindx,I);
   
          % update S and derivative
          S_old=S(cindx,I);
          eval(S_mean_cmd_seq);  S(cindx,I)=f;
          dS(cindx,I)=S(cindx,I)-S_old;        
 
        end   
 
        dSdSN=sum(dS.*dS,1);
        I=find(dSdSN > S_tolN);         
        dSdS=sum(dSdSN);

      end      

  end % switch solver

  % calculate G
  switch solver 
    case 'beliefprop2' 
      for i=1:N
        diagOmega=diag(Omega(:,i)); 
        G(:,:,i)=eye(M)+diagOmega*squeeze(G(:,:,i));
      end
    otherwise
      I=1:N; hcav=J*S-V.*S; eval(S_mean_cmd_all); Omega=dfdg./(1+dfdg.*V);
      for i=1:N
        diagOmega=diag(Omega(:,i)); 
        Gi=inv(eye(M)-diagOmega*J); diagG(:,i)=diag(Gi); G(:,:,i)=Gi;  
      end
  end  

  % calculate chi
  for i=1:N
    diagOmega=diag(Omega(:,i));
    chi(:,:,i)=squeeze(G(:,:,i))*diagOmega;
  end 

  if draw
    hcav=J*S-V.*S; I=1:N; eval(S_mean_cmd_all);     
    for i=1:N
      diagchi(:,i)=diag(squeeze(chi(:,:,i)));
    end    
    zerochi=(dfdg<eps & diagchi<eps);
    dVdV=sum(sum(((1-zerochi).*(1./(dfdg+zerochi)-1./(diagchi+zerochi))).^2));
    fprintf(' Iteration %d\n E-step: S_ite = %d   dSdS = %g   dVdV = %g\n',EM_step,S_ite,dSdS,dVdV);    
    if draw==2
      loglikelihood=loglikelihood_f(X,f,dfdg,h,hcav,V_0,V,J,G,likelihood_cmd,...
                                    likelihood_arg,Sigma_prior,Sigma,InvSigma);
      fprintf(' loglikelihood = %6.3f\n',loglikelihood);
    end
  end

  find(dSdSN < S_tol);
  if ~isempty(I) 
    NI=length(I);
    for i=1:M  % set up matrices for M-step.
      for j=i:M
        tracechi(i,j)=sum(squeeze(chi(i,j,I)));
        tracechi(j,i)=tracechi(i,j);
        traceSS(i,j)=sum(squeeze(chi(i,j,I))+S(i,I)'.*S(j,I)');
        traceSS(j,i)=traceSS(i,j);
      end;
    end;
    eval(A_cmd);     % e.g. A=A_positive(traceSS,X(:,I)*(S(:,I))',A);   this is correct EM as long as A is independent of Sigma.
    eval(Sigma_cmd);
  else  
    disp('No samples converged - exiting'); exitflag=0; break
  end 

  % termination criteria 
  A_norm=norm(abs(A),1); dA_rel=abs(A_norm-A_norm_old)/A_norm; A_norm_old=A_norm;
  Sigma_norm=norm(abs(Sigma),1); dSigma_rel=abs(Sigma_norm-Sigma_norm_old)/Sigma_norm; Sigma_norm_old=Sigma_norm;

  if draw
    fprintf(' M-step: d|A|/|A| = %g   d|Sigma|/|Sigma| = %g\n',dA_rel,dSigma_rel);
    if dim_Sigma==D 
      fprintf(' noise variance = %g \n\n',sum(Sigma)/D);
    else       
      fprintf(' noise variance = %g \n\n',trace(Sigma)/size(Sigma,1));
    end     
  end

  % update parameters of the prior 

end % EM_step

% calculate loglikelihood
hcav=J*S-V.*S; I=1:N; eval(S_mean_cmd_all);
loglikelihood=loglikelihood_f(X,f,dfdg,h,hcav,V_0,V,J,G,likelihood_cmd,...
                              likelihood_arg,Sigma_prior,Sigma,InvSigma); 
if dim_Sigma==D 
  fprintf(' %d EM steps, loglikelihood = %g, noise variance = %g\n',...
          EM_step,loglikelihood,sum(Sigma)/D);   
else       
  fprintf(' %d EM steps, loglikelihood = %g, noise variance = %g\n',...
          EM_step,loglikelihood,trace(Sigma)/size(Sigma,1)); 
end   
loglikelihood=real(loglikelihood); % small imaginary contributions might occur

% return variables sorted according to 'energy'  
[Y,I]=sort((sum(A.^2,1))'.*sum(S.^2,2));
I=flipud(I);
S=S(I,:); A=A(:,I); 
for i=1:N
  chi(:,:,i)=squeeze(chi(I,I,i));
end 
if EM_step==max_ite exitflag=1;
elseif (dA_rel<tol & dSigma_rel<tol) exitflag=2; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        end of main program                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          help functions                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    mixing matrix (A) functions                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = A_free(traceSS,XSt,A) 
A=XSt*inv(traceSS);

function [A] = A_positive(traceSS,XSt,A) 

A=A+10^-3*(A<eps);
amp=XSt./(A*traceSS);
Aneg=any(any(amp<eps));
KT_max_ite=200; KT_ite=0;
while ~Aneg & KT_ite<KT_max_ite 
  KT_ite=KT_ite+1;
  A=A.*amp;
  amp=XSt./(A*traceSS);
  Aneg=any(any(amp<eps));
end
if Aneg % use quadratic programming instead
  M=size(traceSS,1);
  D=size(XSt,1);
  options=optimset('Display','off','TolX',10^5*eps,'TolFun',10^5*eps);
  for i=1:D
    B=quadprog(traceSS,-XSt(i,:)',[],[],[],[],zeros(M,1),[],A(i,:)',options);
    A(i,:)=B';
  end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       mean (S) functions                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,df] = bigauss(gamma,lambda)
sigma2=1; mu=1;
lambda=1./(1+lambda*sigma2);
f=tanh(mu*gamma.*lambda);
df=lambda.*(sigma2+mu^2.*lambda.*(1-f.^2));
f=lambda.*(sigma2*gamma+mu*f);

function [f,df] = binary(gamma,lambda)
f=tanh(gamma);
df=1-f.^2;

function [f,df] = binary_01(gamma,lambda)
tmp=exp(0.5*lambda+gamma);
f=1./(1+exp(0.5*lambda-gamma));
df=f.*(1-f);

function [f,df] = combi(gamma,lambda,prior)
eval(prior);

function [f,df] = exponential(gamma,lambda)
eta=1;
erfclimit=-35;
minlambda=10^-4;
lambda=lambda.*(lambda>minlambda)+minlambda.*(lambda<=minlambda);
xi=(gamma-eta)./sqrt(lambda);
cc=(xi>erfclimit);
xi1=xi.*cc;
epfrac=exp(-(xi1.^2)/2)./(Phi(xi1)*sqrt(2*pi));
f=cc.*(xi1+epfrac)./sqrt(lambda);            % need to go to higher order to get asymptotics right
df=(1./lambda+f.*(xi1./sqrt(lambda)-f)); % need to go to higher order to get asymptotics right -fix at some point!!!

function [f,df] = Gauss(gamma,lambda)
f=gamma./(1+lambda);
df=1./(1+lambda);

function [f,df] = heavy_tail(gamma,lambda)
alpha=1; % if changed change also in logZ0_heavy_tail
f=gamma./lambda-alpha*gamma./(2*alpha*lambda+gamma.^2);
df=1./lambda+alpha*(gamma.^2-2*alpha*lambda)./(2*alpha*lambda+gamma.^2).^2;

function [f,df] = heavy_tail_plus_delta(gamma,lambda)
alpha=1; % if changed change also in logZ0_heavy_tail_plus_delta
beta=0.3; % proporation delta if changed change also in logZ0_heavy_tail_plus_delta
Z0ht=exp(0.5*gamma.^2./lambda).*(1+gamma.^2./(2*alpha*lambda)).^(-0.5*alpha);
f=(1-beta)*(gamma./lambda-alpha*gamma./(2*alpha*lambda+gamma.^2))./...
  (beta./Z0ht+(1-beta));
df=(1-beta)./(beta./Z0ht+(1-beta)).*...
   (1./lambda+alpha*(gamma.^2-2*alpha*lambda)./(2*alpha*lambda+gamma.^2).^2+...
   (gamma./lambda-alpha*gamma./(2*alpha*lambda+gamma.^2)).^2)-f.^2;

function [f,df] = Laplace(gamma,lambda)
erfclimit=-25;
eta=1;
minlambda=10^-4;
lambda=lambda.*(lambda>minlambda)+minlambda.*(lambda<=minlambda);
xip=(gamma-eta)./sqrt(lambda);
ccp=(xip>erfclimit);ccpc=not(ccp);
xip1=ccp.*xip;
xim=-(gamma+eta)./sqrt(lambda);
ccm=(xim>erfclimit);ccmc=not(ccm);
xim1=ccm.*xim;
Dp=exp(-(xip1.^2)/2)/sqrt(2*pi);
Dm=exp(-(xim1.^2)/2)/sqrt(2*pi);
kp=Phi(xip1).*Dm;  
km=Phi(xim1).*Dp; 
f=ccp.*ccm.*(xip.*kp-xim.*km)./(sqrt(lambda).*(kp+km))+(-ccpc.*xim+ccmc.*xip)./sqrt(lambda);
df=(ccp.*ccm.*(1+xim.*xip+Dp.*Dm.*(xip+xim)./(kp+km)+sqrt(lambda).*(xip.*km-xim.*kp)./(kp+km).*f)+ccpc+ccmc)./lambda;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       logZ0 functions                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [logZ0] = logZ0_bigauss(gamma,lambda)
sigma2=1; mu=1;
lambda=1./(1+lambda*sigma2);
logZ0terms=0.5*(log(lambda)-mu^2/sigma2+lambda.*(mu^2/sigma2+gamma.^2*sigma2))...
           -log(2)+abs(gamma*mu.*lambda)+log(1+exp(-2*abs(gamma*mu.*lambda))); % log(cosh(gamma*mu.*lambda))
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_binary(gamma,lambda)
logZ0terms=-0.5*lambda-log(2)+abs(gamma)+log(1+exp(-2*abs(gamma)));
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_binary_01(gamma,lambda)
gamma=0.5*(gamma-0.5*lambda);
logZ0terms=gamma-log(2)+abs(gamma)+log(1+exp(-2*abs(gamma)));
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_combi(gamma,lambda,prior)
logZ0=0;
eval(prior);

function [logZ0] = logZ0_exponential(gamma,lambda)
eta=1;
erfclimit=-35;
minlambda=10^-4;
lambda=lambda.*(lambda>minlambda)+minlambda.*(lambda<=minlambda);
xi=(gamma-eta)./sqrt(lambda);
cc=(xi>erfclimit);
xi1=xi.*cc;
logZ0terms=cc.*(log(Phi(xi1))+0.5*log(2*pi)+0.5*xi1.^2)-(1-cc).*log(abs(xi)+cc)+log(eta)-0.5*log(lambda);
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_Gauss(gamma,lambda)
logZ0terms=0.5*gamma.^2./(1+lambda)-0.5*log(1+lambda);
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_heavy_tail(gamma,lambda)
alpha=1; % if changed change also in heavy_tail.m
logZ0terms=0.5*gamma.^2./lambda-0.5*alpha*log(1+gamma.^2./(2*lambda*alpha));
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_heavy_tail_plus_delta(gamma,lambda)
alpha=1; % if changed change also in heavy_tail_plus_delta
beta=0.3; % proporation delta if changed change also in heavy_tail_plus_delta
Z0ht=exp(0.5*gamma.^2./lambda).*(1+gamma.^2./(2*alpha*lambda)).^(-0.5*alpha);
logZ0terms=log(beta+(1-beta)*Z0ht);
logZ0=sum(sum(logZ0terms));

function [logZ0] = logZ0_Laplace(gamma,lambda)
erfclimit=-25;
eta=1;
minlambda=10^-4;
lambda=lambda.*(lambda>minlambda)+minlambda.*(lambda<=minlambda);
xip=(gamma-eta)./sqrt(lambda);
ccp=(xip>erfclimit);ccpc=not(ccp);
xip1=ccp.*xip;
xim=-(gamma+eta)./sqrt(lambda);
ccm=(xim>erfclimit);ccmc=not(ccm);
xim1=ccm.*xim;
Dp=exp(-(xip1.^2)/2)/sqrt(2*pi);
Dm=exp(-(xim1.^2)/2)/sqrt(2*pi);
Phip=Phi(xip1);  
Phim=Phi(xim1); 
logZ0terms=log(0.5*eta)+0.5*log(2*pi./lambda)+... 
ccp.*ccm.*(0.5*xip1.^2+0.5*xim1.^2+log(Dm.*Phip+Dp.*Phim)-0.5*log(2*pi))+...
ccp.*ccmc.*(0.5*xip1.^2+log(Phip+Dp./(abs(xim)+ccm)))+...
ccpc.*ccm.*(0.5*xim1.^2+log(Phim+Dm./(abs(xip)+ccp)))+...
ccpc.*ccmc.*(log(1./(abs(xip)+ccp)+1./(abs(xim)+ccm))-0.5*log(2*pi));
logZ0=sum(sum(logZ0terms));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Phi function                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phi function

function y=Phi(x)
%
z=abs(x/sqrt(2));
t=1.0./(1.0+0.5*z);
y=0.5*t.*exp(-z.*z-1.26551223+t.*(1.00002368+...
    t.*(0.37409196+t.*(0.09678418+t.*(-0.18628806+t.*(0.27886807+...
    t.*(-1.13520398+t.*(1.48851587+t.*(-0.82215223+t.*0.17087277)))))))));
y=(x<=0.0).*y+(x>0).*(1.0-y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   loglikelihood function                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loglikelihood]=loglikelihood_f(X,S,dfdg,h,hcav,V_0,V,J,G,likelihood_cmd,likelihood_arg,Sigma_prior,Sigma,InvSigma);
% calculates loglikelihood per S variable
%
[M N]=size(S);
D=size(X,1);
eval(likelihood_cmd);
free=free+sum(sum(hcav.*S))+0.5*sum(sum(V.*S.*S))+0.5*sum(sum(log(1+dfdg.*V)));
for i=1:N
  free=free-0.5*S(:,i)'*J*S(:,i)-0.5*log(det(squeeze(G(:,:,i))));
end
switch size(Sigma,1)*size(Sigma,2)
  case 1   % 'isotropic'
    logdetSigma=D*log(Sigma);  
    XInvSigmaX=InvSigma*sum(sum(X.^2));
  case D   % 'diagonal'
    logdetSigma=sum(log(Sigma)); 
    XInvSigmaX=InvSigma'*sum(X.^2,2);
  case D^2 % 'free'
    [L,U]=lu(Sigma);
    logdetSigma=sum(log(abs(diag(U)))); % assuming positive determinant
    XInvSigmaX=sum(sum(X.*(InvSigma*X)));
end
free=free+0.5*N*(D*log(2*pi)+logdetSigma)+0.5*XInvSigmaX;
loglikelihood=-free/N;

% initialization of method, A, S and Sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initialize method                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prior,par]=method_init(prior,par)
try 
  switch prior.method
    case 'constant'
      prior.A='constant'; prior.Sigma='constant';
    case 'fa'
      prior.A='free'; prior.S='Gauss'; prior.Sigma='diagonal';
    case 'neg_kurtosis'
      prior.S='bigauss';
    case 'pos_kurtosis'
      prior.S='heavy_tail'; 
    case 'positive'
      prior.S='exponential'; prior.A='positive';
    case 'ppca' % probabilistic PCA
      prior.A='free'; prior.S='Gauss'; prior.Sigma='isotropic';
  end
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initialize A                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,A_cmd,A_prior,M]=A_init(prior,par,D);
try A_prior=prior.A; catch A_prior='free'; end                % default mixing matrix
switch A_prior
  case 'constant'
    A_cmd='';  
  case 'MacKay'
    A_cmd=sprintf('A=A_free(traceSS-tracechi,X(:,I)*(S(:,I))''-A*tracechi,A);');
  otherwise
    A_cmd=sprintf('A=A_%s(traceSS,X(:,I)*(S(:,I))'',A);',A_prior);
end
try
  A=par.A_init;
  M=size(A,2);
catch
  try M=par.sources; catch M=D; end                           % quadratic mixing matrix is default 
  aMD=abs(D-M); scale=1; A=scale*randn(D,M)/10;                                        
  if (D>M)                                                    % default is asymmetric toeplitz plus Gaussian noise 
    a=(M+1:-1:M-D+2)/(M+1); A=A+scale*toeplitz(a.*(a>0),(1:M)<=1);
  else 
    a=(D+1:-1:D-M+2)/(D+1); A=A+scale*toeplitz((1:D)<=1,a.*(a>0));
  end
  if strcmp('positive',A_prior) A=A.*sign(A);
  else A=sign(rand(D,M)-0.5).*A; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initialize S                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,S_mean_cmd_all,S_mean_cmd_seq,likelihood_cmd,...
          likelihood_arg,S_arg_seq,cS_prior,cM,S_prior]=S_init(prior,par,M,N);
try S_prior=prior.S; catch S_prior='heavy_tail'; end          % default prior
switch S_prior                                                % set up combi prior
  case 'combi'
    cS_prior={};
    for i=1:50 % maximum different priors
      try S_arg=sprintf('cS_prior{i}=prior.S%i; cM(i)=prior.M%i;',i,i); eval(S_arg); catch break; end 
    end
    if i==1 
      S_prior='heavy_tail';  
      S_mean_cmd_seq=sprintf('[f,dfdg] = %s(h(cindx,I)+hcav,V_0(cindx,I)-V(cindx,I));',S_prior);
      S_mean_cmd_all=sprintf('[f,dfdg] = %s(h(:,I)+hcav,V_0(:,I)-V(:,I));',S_prior); 
      likelihood_cmd=sprintf('free = -logZ0_%s(h+hcav,V_0-V);',S_prior);
    end
    bindx=1; S_arg=''; S_arg_seq={}; likelihood_arg='';
    for j=1:i-1
     eindx=min(M,bindx+cM(j)-1);
     if bindx>M break; end
     if (j==i-1 | bindx+cM(j)>M) eindx=M; end
     S_arg=sprintf('%s[f(%i:%i,:),df(%i:%i,:)] = %s(gamma(%i:%i,:),lambda(%i:%i,:));\n',...
                   S_arg,bindx,eindx,bindx,eindx,cS_prior{j},bindx,eindx,bindx,eindx);
     likelihood_arg=sprintf('%slogZ0 = logZ0 + logZ0_%s(gamma(%i:%i,:),lambda(%i:%i,:));\n',...
                            likelihood_arg,cS_prior{j},bindx,eindx,bindx,eindx);
     for k=bindx:eindx S_arg_seq{k}=sprintf('[f,df] = %s(gamma,lambda);',cS_prior{j}); end
     bindx=bindx+cM(j); 
    end  
    S_mean_cmd_seq=sprintf('[f,dfdg] = combi(h(cindx,I)+hcav,V_0(cindx,I)-V(cindx,I),S_arg_seq{cindx});');
    S_mean_cmd_all=sprintf('[f,dfdg] = combi(h(:,I)+hcav,V_0(:,I)-V(:,I),S_arg);'); 
    likelihood_cmd=sprintf('free = -logZ0_combi(h+hcav,V_0-V,likelihood_arg);');
  otherwise  
    likelihood_arg=''; cS_prior=''; cM=0; S_arg_seq='';
    S_mean_cmd_seq=sprintf('[f,dfdg] = %s(h(cindx,I)+hcav,V_0(cindx,I)-V(cindx,I));',S_prior);
    S_mean_cmd_all=sprintf('[f,dfdg] = %s(h(:,I)+hcav,V_0(:,I)-V(:,I));',S_prior); 
    likelihood_cmd=sprintf('free = -logZ0_%s(h+hcav,V_0-V);',S_prior);
end
try
  S=par.S_init;
catch
  S=zeros(M,N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initialize Sigma                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sigma,Sigma_cmd,Sigma_eta,Sigma_prior,Sigma_rel]=Sigma_init(prior,par,A_prior,S_prior,X,A,S)  
try Sigma_eta=par.Sigma_eta; catch Sigma_eta=1; end           % learning rate for Sigma
try 
  Sigma_prior=prior.Sigma;
catch
  Sigma_prior='isotropic';                                      % default prior
end
switch Sigma_prior
  case 'constant' 
    Sigma_cmd='';
  case 'isotropic'
    Sigma_cmd=sprintf('Sigma = (1-Sigma_eta)*Sigma+Sigma_eta*(sum(sum((X(:,I)-A*S(:,I)).^2))+sum(sum((A*tracechi).*A)))/(NI*D);');    
  case 'diagonal'
    Sigma_cmd=sprintf('Sigma = (1-Sigma_eta)*Sigma+Sigma_eta*(sum((X(:,I)-A*S(:,I)).^2,2)+sum((A*tracechi).*A,2))/NI;');    
  case 'free'
    Sigma_cmd=sprintf('Sigma = (1-Sigma_eta)*Sigma+Sigma_eta*((X(:,I)-A*S(:,I))*(X(:,I)-A*S(:,I))''+A*tracechi*A'')/NI;');      
end
if strcmp('free',Sigma_prior) & strcmp('positive',A_prior)
  fprintf('   Prior analyse: uses isotropic noise approximation in A-update.\n   Full update is computationally expensive. Not implemented yet.');   
end
try 
  Sigma_rel=par.Sigma_rel;
catch
  if strcmp(S_prior,'heavy_tail') % | strcmp(S_prior,'heavy_tail_plus_delta')
    Sigma_rel=1;                                            % default relative noise level 
  else
    Sigma_rel=1;
  end
end 
[D,N]=size(X);
try
  Sigma=par.Sigma_init;  
catch                               
  switch Sigma_prior
    case {'isotropic','constant'}
      Sigma=Sigma_rel*sum(sum((X-A*S).^2))/(D*N);   
    case 'diagonal'
      Sigma=Sigma_rel*sum((X-A*S).^2,2)/N; 
    case 'free'
      Sigma=Sigma_rel*(X-A*S)*(X-A*S)'/N; 
  end
end

