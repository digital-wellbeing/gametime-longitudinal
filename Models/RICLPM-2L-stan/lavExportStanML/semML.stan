functions{

  vector[] sem_mean(vector[] alpha, real[,,] B, real[,,] gamma, int[] g, int k, int Ng, int gamind, real[,] meanx){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(gamind == 1){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (alpha[j] + to_matrix(gamma[,,j]) * to_vector(meanx[,j]));

      } else {
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * alpha[j];
      }
    }

    return evlv;
  }

  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;
    nidx = nov + nlvno0;

    iden = diag_matrix(rep_vector(1.0, k));

    if(nlvno0 > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:nidx){
          psivecinv[j,idx[i]] = 1/psi[idx[i],idx[i],j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,idx,idx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	  ldetcomp[j] = -2 * ldetcomp[j] + sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
	psimat[j] = to_matrix(psi[,,j]) + to_matrix(psi[,,j])' - diag_matrix(diagonal(to_matrix(psi[,,j])));

	ldetcomp[j] = log_determinant(psimat[j,idx,idx]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] - 2 * log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	}

	psimatinv[j] = psimat[j];
	psimatinv[j,1:nidx,1:nidx] = inverse_spd(psimat[j,idx,idx]);
        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,1:nidx,1:nidx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec[idx] - evlv[g[i],idx])' * siginv[g[i],1:nidx,1:nidx] * (xvec[idx] - evlv[g[i],idx]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }

  matrix fill_lower(matrix x){
    matrix[rows(x),cols(x)] newx;

    newx = x;
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        newx[j,i] = x[i,j];
      }
    }
    return newx;
  }
}

data{
  int N;
  int g[N];
  int lvind[8];
  int etaind[8];
  real sampmean[6,10];
  real meanx[6,10];
  int dummyov[2];
  int dummylv[2];
  vector[6] y[N];
  real lambdaframe[6,8,10];
  real thetaframe[6,6,10];
  real psiframe[8,8,10];
  real betaframe[8,8,10];
  real nuframe[6,1,10];
  real alphaframe[8,1,10];
  real lvrhoframe[8,8,10];
  int<lower=1> N_1;  // number of participants per game
  int<lower=1> N_2;  // number of games
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
}

parameters{
  vector<lower=0>[8] psifree;
  vector[4] betafree;
  vector[6] nufree;
  vector<lower=0,upper=1>[4] lvrhofree;
  matrix[N, 8] etavec;
  vector<lower=0>[M_1] sd_1;  // game-level standard deviations
  vector[N_2] z_1[M_1];  // standardized game-level effects
}

transformed parameters{
  real lambda[6,8,10];
  real theta[6,6,10];
  matrix[6,6] thetld[10];
  real psi[8,8,10];
  real beta[8,8,10];
  real nu[6,1,10];
  real alpha[8,1,10];
  real lvrho[8,8,10];
  real mu[N,6];
  matrix[N,8] eta;

  eta = rep_matrix(0, N, 8);

  lambda = lambdaframe;
  theta = thetaframe;
  psi = psiframe;
  beta = betaframe;
  nu = nuframe;
  alpha = alphaframe;
  lvrho = lvrhoframe;

  // Alot of unessecary repetition here
  lambda[1,1,1] = 1;
  lambda[2,1,1] = 1;
  lambda[3,1,1] = 1;
  lambda[4,2,1] = 1;
  lambda[5,2,1] = 1;
  lambda[6,2,1] = 1;
  lambda[1,3,1] = 1;
  lambda[2,4,1] = 1;
  lambda[3,5,1] = 1;
  lambda[4,6,1] = 1;
  lambda[5,7,1] = 1;
  lambda[6,8,1] = 1;
  psi[1,1,1] = pow(psifree[1],2);
  psi[2,2,1] = pow(psifree[2],2);
  psi[3,3,1] = pow(psifree[3],2);
  psi[6,6,1] = pow(psifree[4],2);
  psi[4,4,1] = pow(psifree[5],2);
  psi[7,7,1] = pow(psifree[6],2);
  psi[5,5,1] = pow(psifree[7],2);
  psi[8,8,1] = pow(psifree[8],2);
  theta[1,1,1] = 0.001;
  theta[2,2,1] = 0.001;
  theta[3,3,1] = 0.001;
  theta[4,4,1] = 0.001;
  theta[5,5,1] = 0.001;
  theta[6,6,1] = 0.001;
  nu[1,1,1] = nufree[1];
  nu[2,1,1] = nufree[2];
  nu[3,1,1] = nufree[3];
  nu[4,1,1] = nufree[4];
  nu[5,1,1] = nufree[5];
  nu[6,1,1] = nufree[6];
  alpha[1,1,1] = 0;
  alpha[2,1,1] = 0;
  alpha[3,1,1] = 0;
  alpha[4,1,1] = 0;
  alpha[5,1,1] = 0;
  alpha[6,1,1] = 0;
  alpha[7,1,1] = 0;
  alpha[8,1,1] = 0;
  lvrho[3,6,1] = -1 + 2*lvrhofree[1];
  lvrho[4,7,1] = -1 + 2*lvrhofree[2];
  lvrho[5,8,1] = -1 + 2*lvrhofree[3];
  lvrho[1,2,1] = -1 + 2*lvrhofree[4];
  psi[3,6,1] = lvrho[3,6,1] * sqrt(psi[3,3,1] * psi[6,6,1]);
  psi[4,7,1] = lvrho[4,7,1] * sqrt(psi[4,4,1] * psi[7,7,1]);
  psi[5,8,1] = lvrho[5,8,1] * sqrt(psi[5,5,1] * psi[8,8,1]);
  psi[1,2,1] = lvrho[1,2,1] * sqrt(psi[1,1,1] * psi[2,2,1]);
  lambda[1,1,2] = 1;
  lambda[2,1,2] = 1;
  lambda[3,1,2] = 1;
  lambda[4,2,2] = 1;
  lambda[5,2,2] = 1;
  lambda[6,2,2] = 1;
  lambda[1,3,2] = 1;
  lambda[2,4,2] = 1;
  lambda[3,5,2] = 1;
  lambda[4,6,2] = 1;
  lambda[5,7,2] = 1;
  lambda[6,8,2] = 1;
  psi[1,1,2] = pow(psifree[1],2);
  psi[2,2,2] = pow(psifree[2],2);
  psi[3,3,2] = pow(psifree[3],2);
  psi[6,6,2] = pow(psifree[4],2);
  psi[4,4,2] = pow(psifree[5],2);
  psi[7,7,2] = pow(psifree[6],2);
  psi[5,5,2] = pow(psifree[7],2);
  psi[8,8,2] = pow(psifree[8],2);
  theta[1,1,2] = 0.001;
  theta[2,2,2] = 0.001;
  theta[3,3,2] = 0.001;
  theta[4,4,2] = 0.001;
  theta[5,5,2] = 0.001;
  theta[6,6,2] = 0.001;
  nu[1,1,2] = nufree[1];
  nu[2,1,2] = nufree[2];
  nu[3,1,2] = nufree[3];
  nu[4,1,2] = nufree[4];
  nu[5,1,2] = nufree[5];
  nu[6,1,2] = nufree[6];
  alpha[1,1,2] = 0;
  alpha[2,1,2] = 0;
  alpha[3,1,2] = 0;
  alpha[4,1,2] = 0;
  alpha[5,1,2] = 0;
  alpha[6,1,2] = 0;
  alpha[7,1,2] = 0;
  alpha[8,1,2] = 0;
  lvrho[3,6,2] = -1 + 2*lvrhofree[1];
  lvrho[4,7,2] = -1 + 2*lvrhofree[2];
  lvrho[5,8,2] = -1 + 2*lvrhofree[3];
  lvrho[1,2,2] = -1 + 2*lvrhofree[4];
  psi[3,6,2] = lvrho[3,6,2] * sqrt(psi[3,3,2] * psi[6,6,2]);
  psi[4,7,2] = lvrho[4,7,2] * sqrt(psi[4,4,2] * psi[7,7,2]);
  psi[5,8,2] = lvrho[5,8,2] * sqrt(psi[5,5,2] * psi[8,8,2]);
  psi[1,2,2] = lvrho[1,2,2] * sqrt(psi[1,1,2] * psi[2,2,2]);
  lambda[1,1,3] = 1;
  lambda[2,1,3] = 1;
  lambda[3,1,3] = 1;
  lambda[4,2,3] = 1;
  lambda[5,2,3] = 1;
  lambda[6,2,3] = 1;
  lambda[1,3,3] = 1;
  lambda[2,4,3] = 1;
  lambda[3,5,3] = 1;
  lambda[4,6,3] = 1;
  lambda[5,7,3] = 1;
  lambda[6,8,3] = 1;
  psi[1,1,3] = pow(psifree[1],2);
  psi[2,2,3] = pow(psifree[2],2);
  psi[3,3,3] = pow(psifree[3],2);
  psi[6,6,3] = pow(psifree[4],2);
  psi[4,4,3] = pow(psifree[5],2);
  psi[7,7,3] = pow(psifree[6],2);
  psi[5,5,3] = pow(psifree[7],2);
  psi[8,8,3] = pow(psifree[8],2);
  theta[1,1,3] = 0.001;
  theta[2,2,3] = 0.001;
  theta[3,3,3] = 0.001;
  theta[4,4,3] = 0.001;
  theta[5,5,3] = 0.001;
  theta[6,6,3] = 0.001;
  nu[1,1,3] = nufree[1];
  nu[2,1,3] = nufree[2];
  nu[3,1,3] = nufree[3];
  nu[4,1,3] = nufree[4];
  nu[5,1,3] = nufree[5];
  nu[6,1,3] = nufree[6];
  alpha[1,1,3] = 0;
  alpha[2,1,3] = 0;
  alpha[3,1,3] = 0;
  alpha[4,1,3] = 0;
  alpha[5,1,3] = 0;
  alpha[6,1,3] = 0;
  alpha[7,1,3] = 0;
  alpha[8,1,3] = 0;
  lvrho[3,6,3] = -1 + 2*lvrhofree[1];
  lvrho[4,7,3] = -1 + 2*lvrhofree[2];
  lvrho[5,8,3] = -1 + 2*lvrhofree[3];
  lvrho[1,2,3] = -1 + 2*lvrhofree[4];
  psi[3,6,3] = lvrho[3,6,3] * sqrt(psi[3,3,3] * psi[6,6,3]);
  psi[4,7,3] = lvrho[4,7,3] * sqrt(psi[4,4,3] * psi[7,7,3]);
  psi[5,8,3] = lvrho[5,8,3] * sqrt(psi[5,5,3] * psi[8,8,3]);
  psi[1,2,3] = lvrho[1,2,3] * sqrt(psi[1,1,3] * psi[2,2,3]);
  lambda[1,1,4] = 1;
  lambda[2,1,4] = 1;
  lambda[3,1,4] = 1;
  lambda[4,2,4] = 1;
  lambda[5,2,4] = 1;
  lambda[6,2,4] = 1;
  lambda[1,3,4] = 1;
  lambda[2,4,4] = 1;
  lambda[3,5,4] = 1;
  lambda[4,6,4] = 1;
  lambda[5,7,4] = 1;
  lambda[6,8,4] = 1;
  psi[1,1,4] = pow(psifree[1],2);
  psi[2,2,4] = pow(psifree[2],2);
  psi[3,3,4] = pow(psifree[3],2);
  psi[6,6,4] = pow(psifree[4],2);
  psi[4,4,4] = pow(psifree[5],2);
  psi[7,7,4] = pow(psifree[6],2);
  psi[5,5,4] = pow(psifree[7],2);
  psi[8,8,4] = pow(psifree[8],2);
  theta[1,1,4] = 0.001;
  theta[2,2,4] = 0.001;
  theta[3,3,4] = 0.001;
  theta[4,4,4] = 0.001;
  theta[5,5,4] = 0.001;
  theta[6,6,4] = 0.001;
  nu[1,1,4] = nufree[1];
  nu[2,1,4] = nufree[2];
  nu[3,1,4] = nufree[3];
  nu[4,1,4] = nufree[4];
  nu[5,1,4] = nufree[5];
  nu[6,1,4] = nufree[6];
  alpha[1,1,4] = 0;
  alpha[2,1,4] = 0;
  alpha[3,1,4] = 0;
  alpha[4,1,4] = 0;
  alpha[5,1,4] = 0;
  alpha[6,1,4] = 0;
  alpha[7,1,4] = 0;
  alpha[8,1,4] = 0;
  lvrho[3,6,4] = -1 + 2*lvrhofree[1];
  lvrho[4,7,4] = -1 + 2*lvrhofree[2];
  lvrho[5,8,4] = -1 + 2*lvrhofree[3];
  lvrho[1,2,4] = -1 + 2*lvrhofree[4];
  psi[3,6,4] = lvrho[3,6,4] * sqrt(psi[3,3,4] * psi[6,6,4]);
  psi[4,7,4] = lvrho[4,7,4] * sqrt(psi[4,4,4] * psi[7,7,4]);
  psi[5,8,4] = lvrho[5,8,4] * sqrt(psi[5,5,4] * psi[8,8,4]);
  psi[1,2,4] = lvrho[1,2,4] * sqrt(psi[1,1,4] * psi[2,2,4]);
  lambda[1,1,5] = 1;
  lambda[2,1,5] = 1;
  lambda[3,1,5] = 1;
  lambda[4,2,5] = 1;
  lambda[5,2,5] = 1;
  lambda[6,2,5] = 1;
  lambda[1,3,5] = 1;
  lambda[2,4,5] = 1;
  lambda[3,5,5] = 1;
  lambda[4,6,5] = 1;
  lambda[5,7,5] = 1;
  lambda[6,8,5] = 1;
  psi[1,1,5] = pow(psifree[1],2);
  psi[2,2,5] = pow(psifree[2],2);
  psi[3,3,5] = pow(psifree[3],2);
  psi[6,6,5] = pow(psifree[4],2);
  psi[4,4,5] = pow(psifree[5],2);
  psi[7,7,5] = pow(psifree[6],2);
  psi[5,5,5] = pow(psifree[7],2);
  psi[8,8,5] = pow(psifree[8],2);
  theta[1,1,5] = 0.001;
  theta[2,2,5] = 0.001;
  theta[3,3,5] = 0.001;
  theta[4,4,5] = 0.001;
  theta[5,5,5] = 0.001;
  theta[6,6,5] = 0.001;
  nu[1,1,5] = nufree[1];
  nu[2,1,5] = nufree[2];
  nu[3,1,5] = nufree[3];
  nu[4,1,5] = nufree[4];
  nu[5,1,5] = nufree[5];
  nu[6,1,5] = nufree[6];
  alpha[1,1,5] = 0;
  alpha[2,1,5] = 0;
  alpha[3,1,5] = 0;
  alpha[4,1,5] = 0;
  alpha[5,1,5] = 0;
  alpha[6,1,5] = 0;
  alpha[7,1,5] = 0;
  alpha[8,1,5] = 0;
  lvrho[3,6,5] = -1 + 2*lvrhofree[1];
  lvrho[4,7,5] = -1 + 2*lvrhofree[2];
  lvrho[5,8,5] = -1 + 2*lvrhofree[3];
  lvrho[1,2,5] = -1 + 2*lvrhofree[4];
  psi[3,6,5] = lvrho[3,6,5] * sqrt(psi[3,3,5] * psi[6,6,5]);
  psi[4,7,5] = lvrho[4,7,5] * sqrt(psi[4,4,5] * psi[7,7,5]);
  psi[5,8,5] = lvrho[5,8,5] * sqrt(psi[5,5,5] * psi[8,8,5]);
  psi[1,2,5] = lvrho[1,2,5] * sqrt(psi[1,1,5] * psi[2,2,5]);
  lambda[1,1,6] = 1;
  lambda[2,1,6] = 1;
  lambda[3,1,6] = 1;
  lambda[4,2,6] = 1;
  lambda[5,2,6] = 1;
  lambda[6,2,6] = 1;
  lambda[1,3,6] = 1;
  lambda[2,4,6] = 1;
  lambda[3,5,6] = 1;
  lambda[4,6,6] = 1;
  lambda[5,7,6] = 1;
  lambda[6,8,6] = 1;
  psi[1,1,6] = pow(psifree[1],2);
  psi[2,2,6] = pow(psifree[2],2);
  psi[3,3,6] = pow(psifree[3],2);
  psi[6,6,6] = pow(psifree[4],2);
  psi[4,4,6] = pow(psifree[5],2);
  psi[7,7,6] = pow(psifree[6],2);
  psi[5,5,6] = pow(psifree[7],2);
  psi[8,8,6] = pow(psifree[8],2);
  theta[1,1,6] = 0.001;
  theta[2,2,6] = 0.001;
  theta[3,3,6] = 0.001;
  theta[4,4,6] = 0.001;
  theta[5,5,6] = 0.001;
  theta[6,6,6] = 0.001;
  nu[1,1,6] = nufree[1];
  nu[2,1,6] = nufree[2];
  nu[3,1,6] = nufree[3];
  nu[4,1,6] = nufree[4];
  nu[5,1,6] = nufree[5];
  nu[6,1,6] = nufree[6];
  alpha[1,1,6] = 0;
  alpha[2,1,6] = 0;
  alpha[3,1,6] = 0;
  alpha[4,1,6] = 0;
  alpha[5,1,6] = 0;
  alpha[6,1,6] = 0;
  alpha[7,1,6] = 0;
  alpha[8,1,6] = 0;
  lvrho[3,6,6] = -1 + 2*lvrhofree[1];
  lvrho[4,7,6] = -1 + 2*lvrhofree[2];
  lvrho[5,8,6] = -1 + 2*lvrhofree[3];
  lvrho[1,2,6] = -1 + 2*lvrhofree[4];
  psi[3,6,6] = lvrho[3,6,6] * sqrt(psi[3,3,6] * psi[6,6,6]);
  psi[4,7,6] = lvrho[4,7,6] * sqrt(psi[4,4,6] * psi[7,7,6]);
  psi[5,8,6] = lvrho[5,8,6] * sqrt(psi[5,5,6] * psi[8,8,6]);
  psi[1,2,6] = lvrho[1,2,6] * sqrt(psi[1,1,6] * psi[2,2,6]);
  lambda[1,1,7] = 1;
  lambda[2,1,7] = 1;
  lambda[3,1,7] = 1;
  lambda[4,2,7] = 1;
  lambda[5,2,7] = 1;
  lambda[6,2,7] = 1;
  lambda[1,3,7] = 1;
  lambda[2,4,7] = 1;
  lambda[3,5,7] = 1;
  lambda[4,6,7] = 1;
  lambda[5,7,7] = 1;
  lambda[6,8,7] = 1;
  psi[1,1,7] = pow(psifree[1],2);
  psi[2,2,7] = pow(psifree[2],2);
  psi[3,3,7] = pow(psifree[3],2);
  psi[6,6,7] = pow(psifree[4],2);
  psi[4,4,7] = pow(psifree[5],2);
  psi[7,7,7] = pow(psifree[6],2);
  psi[5,5,7] = pow(psifree[7],2);
  psi[8,8,7] = pow(psifree[8],2);
  theta[1,1,7] = 0.001;
  theta[2,2,7] = 0.001;
  theta[3,3,7] = 0.001;
  theta[4,4,7] = 0.001;
  theta[5,5,7] = 0.001;
  theta[6,6,7] = 0.001;
  nu[1,1,7] = nufree[1];
  nu[2,1,7] = nufree[2];
  nu[3,1,7] = nufree[3];
  nu[4,1,7] = nufree[4];
  nu[5,1,7] = nufree[5];
  nu[6,1,7] = nufree[6];
  alpha[1,1,7] = 0;
  alpha[2,1,7] = 0;
  alpha[3,1,7] = 0;
  alpha[4,1,7] = 0;
  alpha[5,1,7] = 0;
  alpha[6,1,7] = 0;
  alpha[7,1,7] = 0;
  alpha[8,1,7] = 0;
  lvrho[3,6,7] = -1 + 2*lvrhofree[1];
  lvrho[4,7,7] = -1 + 2*lvrhofree[2];
  lvrho[5,8,7] = -1 + 2*lvrhofree[3];
  lvrho[1,2,7] = -1 + 2*lvrhofree[4];
  psi[3,6,7] = lvrho[3,6,7] * sqrt(psi[3,3,7] * psi[6,6,7]);
  psi[4,7,7] = lvrho[4,7,7] * sqrt(psi[4,4,7] * psi[7,7,7]);
  psi[5,8,7] = lvrho[5,8,7] * sqrt(psi[5,5,7] * psi[8,8,7]);
  psi[1,2,7] = lvrho[1,2,7] * sqrt(psi[1,1,7] * psi[2,2,7]);
  lambda[1,1,8] = 1;
  lambda[2,1,8] = 1;
  lambda[3,1,8] = 1;
  lambda[4,2,8] = 1;
  lambda[5,2,8] = 1;
  lambda[6,2,8] = 1;
  lambda[1,3,8] = 1;
  lambda[2,4,8] = 1;
  lambda[3,5,8] = 1;
  lambda[4,6,8] = 1;
  lambda[5,7,8] = 1;
  lambda[6,8,8] = 1;
  psi[1,1,8] = pow(psifree[1],2);
  psi[2,2,8] = pow(psifree[2],2);
  psi[3,3,8] = pow(psifree[3],2);
  psi[6,6,8] = pow(psifree[4],2);
  psi[4,4,8] = pow(psifree[5],2);
  psi[7,7,8] = pow(psifree[6],2);
  psi[5,5,8] = pow(psifree[7],2);
  psi[8,8,8] = pow(psifree[8],2);
  theta[1,1,8] = 0.001;
  theta[2,2,8] = 0.001;
  theta[3,3,8] = 0.001;
  theta[4,4,8] = 0.001;
  theta[5,5,8] = 0.001;
  theta[6,6,8] = 0.001;
  nu[1,1,8] = nufree[1];
  nu[2,1,8] = nufree[2];
  nu[3,1,8] = nufree[3];
  nu[4,1,8] = nufree[4];
  nu[5,1,8] = nufree[5];
  nu[6,1,8] = nufree[6];
  alpha[1,1,8] = 0;
  alpha[2,1,8] = 0;
  alpha[3,1,8] = 0;
  alpha[4,1,8] = 0;
  alpha[5,1,8] = 0;
  alpha[6,1,8] = 0;
  alpha[7,1,8] = 0;
  alpha[8,1,8] = 0;
  lvrho[3,6,8] = -1 + 2*lvrhofree[1];
  lvrho[4,7,8] = -1 + 2*lvrhofree[2];
  lvrho[5,8,8] = -1 + 2*lvrhofree[3];
  lvrho[1,2,8] = -1 + 2*lvrhofree[4];
  psi[3,6,8] = lvrho[3,6,8] * sqrt(psi[3,3,8] * psi[6,6,8]);
  psi[4,7,8] = lvrho[4,7,8] * sqrt(psi[4,4,8] * psi[7,7,8]);
  psi[5,8,8] = lvrho[5,8,8] * sqrt(psi[5,5,8] * psi[8,8,8]);
  psi[1,2,8] = lvrho[1,2,8] * sqrt(psi[1,1,8] * psi[2,2,8]);
  lambda[1,1,9] = 1;
  lambda[2,1,9] = 1;
  lambda[3,1,9] = 1;
  lambda[4,2,9] = 1;
  lambda[5,2,9] = 1;
  lambda[6,2,9] = 1;
  lambda[1,3,9] = 1;
  lambda[2,4,9] = 1;
  lambda[3,5,9] = 1;
  lambda[4,6,9] = 1;
  lambda[5,7,9] = 1;
  lambda[6,8,9] = 1;
  psi[1,1,9] = pow(psifree[1],2);
  psi[2,2,9] = pow(psifree[2],2);
  psi[3,3,9] = pow(psifree[3],2);
  psi[6,6,9] = pow(psifree[4],2);
  psi[4,4,9] = pow(psifree[5],2);
  psi[7,7,9] = pow(psifree[6],2);
  psi[5,5,9] = pow(psifree[7],2);
  psi[8,8,9] = pow(psifree[8],2);
  theta[1,1,9] = 0.001;
  theta[2,2,9] = 0.001;
  theta[3,3,9] = 0.001;
  theta[4,4,9] = 0.001;
  theta[5,5,9] = 0.001;
  theta[6,6,9] = 0.001;
  nu[1,1,9] = nufree[1];
  nu[2,1,9] = nufree[2];
  nu[3,1,9] = nufree[3];
  nu[4,1,9] = nufree[4];
  nu[5,1,9] = nufree[5];
  nu[6,1,9] = nufree[6];
  alpha[1,1,9] = 0;
  alpha[2,1,9] = 0;
  alpha[3,1,9] = 0;
  alpha[4,1,9] = 0;
  alpha[5,1,9] = 0;
  alpha[6,1,9] = 0;
  alpha[7,1,9] = 0;
  alpha[8,1,9] = 0;
  lvrho[3,6,9] = -1 + 2*lvrhofree[1];
  lvrho[4,7,9] = -1 + 2*lvrhofree[2];
  lvrho[5,8,9] = -1 + 2*lvrhofree[3];
  lvrho[1,2,9] = -1 + 2*lvrhofree[4];
  psi[3,6,9] = lvrho[3,6,9] * sqrt(psi[3,3,9] * psi[6,6,9]);
  psi[4,7,9] = lvrho[4,7,9] * sqrt(psi[4,4,9] * psi[7,7,9]);
  psi[5,8,9] = lvrho[5,8,9] * sqrt(psi[5,5,9] * psi[8,8,9]);
  psi[1,2,9] = lvrho[1,2,9] * sqrt(psi[1,1,9] * psi[2,2,9]);
  lambda[1,1,10] = 1;
  lambda[2,1,10] = 1;
  lambda[3,1,10] = 1;
  lambda[4,2,10] = 1;
  lambda[5,2,10] = 1;
  lambda[6,2,10] = 1;
  lambda[1,3,10] = 1;
  lambda[2,4,10] = 1;
  lambda[3,5,10] = 1;
  lambda[4,6,10] = 1;
  lambda[5,7,10] = 1;
  lambda[6,8,10] = 1;
  psi[1,1,10] = pow(psifree[1],2);
  psi[2,2,10] = pow(psifree[2],2);
  psi[3,3,10] = pow(psifree[3],2);
  psi[6,6,10] = pow(psifree[4],2);
  psi[4,4,10] = pow(psifree[5],2);
  psi[7,7,10] = pow(psifree[6],2);
  psi[5,5,10] = pow(psifree[7],2);
  psi[8,8,10] = pow(psifree[8],2);
  theta[1,1,10] = 0.001;
  theta[2,2,10] = 0.001;
  theta[3,3,10] = 0.001;
  theta[4,4,10] = 0.001;
  theta[5,5,10] = 0.001;
  theta[6,6,10] = 0.001;
  nu[1,1,10] = nufree[1];
  nu[2,1,10] = nufree[2];
  nu[3,1,10] = nufree[3];
  nu[4,1,10] = nufree[4];
  nu[5,1,10] = nufree[5];
  nu[6,1,10] = nufree[6];
  alpha[1,1,10] = 0;
  alpha[2,1,10] = 0;
  alpha[3,1,10] = 0;
  alpha[4,1,10] = 0;
  alpha[5,1,10] = 0;
  alpha[6,1,10] = 0;
  alpha[7,1,10] = 0;
  alpha[8,1,10] = 0;
  lvrho[3,6,10] = -1 + 2*lvrhofree[1];
  lvrho[4,7,10] = -1 + 2*lvrhofree[2];
  lvrho[5,8,10] = -1 + 2*lvrhofree[3];
  lvrho[1,2,10] = -1 + 2*lvrhofree[4];
  psi[3,6,10] = lvrho[3,6,10] * sqrt(psi[3,3,10] * psi[6,6,10]);
  psi[4,7,10] = lvrho[4,7,10] * sqrt(psi[4,4,10] * psi[7,7,10]);
  psi[5,8,10] = lvrho[5,8,10] * sqrt(psi[5,5,10] * psi[8,8,10]);
  psi[1,2,10] = lvrho[1,2,10] * sqrt(psi[1,1,10] * psi[2,2,10]);


  vector[N_2] r_1_1;  // actual group-level effects
  vector[N_2] r_1_2;  // actual group-level effects
  vector[N_2] r_1_3;  // actual group-level effects
  vector[N_2] r_1_4;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_1_2 = (sd_1[2] * (z_1[2]));
  r_1_3 = (sd_1[3] * (z_1[3]));
  r_1_4 = (sd_1[4] * (z_1[4]));

  for(j in 1:N_2) {
    beta[4,3,j] = betafree[1] + r_1_1[j];
    beta[4,6,j] = betafree[2] + r_1_2[j];
    beta[7,3,j] = betafree[3] + r_1_3[j];
    beta[7,6,j] = betafree[4] + r_1_4[j];
    beta[5,4,j] = betafree[1] + r_1_1[j];
    beta[5,7,j] = betafree[2] + r_1_2[j];
    beta[8,4,j] = betafree[3] + r_1_3[j];
    beta[8,7,j] = betafree[4] + r_1_4[j];
  }

  // mu definitions
  for(i in 1:N) {
    eta[i,1:8] = etavec[i];

    mu[i,1] = nu[1,1,g[i]] + lambda[1,1,g[i]]*eta[i,1] + lambda[1,3,g[i]]*eta[i,3];
    mu[i,2] = nu[2,1,g[i]] + lambda[2,1,g[i]]*eta[i,1] + lambda[2,4,g[i]]*eta[i,4];
    mu[i,3] = nu[3,1,g[i]] + lambda[3,1,g[i]]*eta[i,1] + lambda[3,5,g[i]]*eta[i,5];
    mu[i,4] = nu[4,1,g[i]] + lambda[4,2,g[i]]*eta[i,2] + lambda[4,6,g[i]]*eta[i,6];
    mu[i,5] = nu[5,1,g[i]] + lambda[5,2,g[i]]*eta[i,2] + lambda[5,7,g[i]]*eta[i,7];
    mu[i,6] = nu[6,1,g[i]] + lambda[6,2,g[i]]*eta[i,2] + lambda[6,8,g[i]]*eta[i,8];
  }

  for(j in 1:10){
    thetld[j] = fill_lower(to_matrix(theta[,,j]));
    thetld[j] = cholesky_decompose(thetld[j]);
  }

}

model {
  for(i in 1:N) {
    y[i] ~ multi_normal_cholesky(to_vector(mu[i,1:6]), thetld[g[i]]);
  }

  eta ~ sem_lv(alpha, beta, psi, beta, 0, meanx, g, 8, N, 10, 0, 1, 8, etaind, 8);

  // Priors
  target += normal_lpdf(betafree[1] | 0,2);
  target += normal_lpdf(betafree[2] | 0,2);
  target += normal_lpdf(betafree[3] | 0,2);
  target += normal_lpdf(betafree[4] | 0,2);
  target += gamma_lpdf(psifree[1] | 1,.5);
  target += gamma_lpdf(psifree[2] | 1,.5);
  target += gamma_lpdf(psifree[3] | 1,.5);
  target += gamma_lpdf(psifree[4] | 1,.5);
  target += gamma_lpdf(psifree[5] | 1,.5);
  target += gamma_lpdf(psifree[6] | 1,.5);
  target += gamma_lpdf(psifree[7] | 1,.5);
  target += gamma_lpdf(psifree[8] | 1,.5);
  target += normal_lpdf(nufree[1] | 0,32);
  target += normal_lpdf(nufree[2] | 0,32);
  target += normal_lpdf(nufree[3] | 0,32);
  target += normal_lpdf(nufree[4] | 0,32);
  target += normal_lpdf(nufree[5] | 0,32);
  target += normal_lpdf(nufree[6] | 0,32);
  target += beta_lpdf(lvrhofree[1] | 1,1);
  target += beta_lpdf(lvrhofree[2] | 1,1);
  target += beta_lpdf(lvrhofree[3] | 1,1);
  target += beta_lpdf(lvrhofree[4] | 1,1);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_1[2]);
  target += std_normal_lpdf(z_1[3]);
  target += std_normal_lpdf(z_1[4]);
}