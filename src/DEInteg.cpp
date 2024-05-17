//
// Created by miguel on 11/05/2024.
//

#include "../include/DEInteg.h"
#include "../include/sign.h"

/*
%----------------------------------------------------------------------------
%
% Purpose:
%   Numerical integration methods for ordinaray differential equations
%
%   This module provides implemenation of the variable order variable
%   stepsize multistep method of Shampine & Gordon.
%
% Last modified:   2015/08/25   M. Mahooti
%
% Reference:
%
%   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
%   Freeman and Comp., San Francisco (1975).
%
%----------------------------------------------------------------------------
*/
Matrix DEInteg(Matrix (*func)(double, Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y) {
    //% maxnum = 500;
    double twou = 2 * (std::numeric_limits<double>::epsilon());
    double fouru = 4 * (std::numeric_limits<double>::epsilon());

    enum DE_STATE {
        DE_INIT = 1,      // Restart integration
        DE_DONE = 2,      // Successful step
        DE_BADACC = 3,    // Accuracy requirement could not be achieved
        DE_NUMSTEPS = 4,  // Permitted number of steps exceeded
        DE_STIFF = 5,     // Stiff problem suspected
        DE_INVPARAM = 6   // Invalid input parameters
    };

    int State_ = DE_STATE::DE_INIT;
    bool PermitTOUT = true;         //% Allow integration past tout by default
    double told = 0.0;

    //% Powers of two (two(n)=2^n)
    double two[] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};

    double gstr[] = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592,
                     0.00524, 0.00468};

    Matrix yy(n_eqn, 1);    //% Allocate vectors with proper dimension
    Matrix wt(n_eqn, 1);
    Matrix p(n_eqn, 1);
    Matrix yp(n_eqn, 1);
    Matrix phi(n_eqn, 17);
    Matrix g(14, 1);
    Matrix sig(14, 1);
    Matrix rho(14, 1);
    Matrix w(13, 1);
    Matrix alpha(13, 1);
    Matrix beta(13, 1);
    Matrix v(13, 1);
    Matrix psi_(13, 1);

    //% while(true)

    //% Return, if output time equals input time

    if (t == tout) {   //% No integration
        return y;
    }

    //% Test for improper parameters

    double epsilon = fmax(relerr, abserr);

    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 ||
        State_ > DE_STATE::DE_INVPARAM || (State_ != DE_STATE::DE_INIT && t != told)) {
        State_ = DE_STATE::DE_INVPARAM; // Establecer el código de error
        return y; // Salir
    }

    //% On each call set interval of integration and counter for
    //% number of steps. Adjust input error tolerances to define
    //% weight vector for subroutine STEP.

    double del = tout - t;
    double absdel = fabs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT) {
        tend = tout;
    }

    int nostep = 0;
    int kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    double x;
    bool start = false;
    if (State_ == DE_STATE::DE_INIT /*|| !OldPermit || delsgn * del <= 0.0*/) {
        //% On start and restart also set the work variables x and yy(*),
        //% store the direction of integration and initialize the step size
        start = true;
        x = t;
        yy = y;
        double delsgn = sign_(1.0, del);
        double h = sign_(fmax(fouru * fabs(x), fabs(tout - x)), tout - x);
    }

    bool OldPermit;
    while (true) { // % Start step loop

        //   % If already past output point, interpolate solution and return
        if (fabs(x - t) >= absdel) {
            Matrix yout(n_eqn, 0.0);
            Matrix ypout(n_eqn, 0.0);
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            double hi = tout - x;
            int ki = 11 + 1;

            //% Initialize w[*] for computing g[*]
            for (int i = 1; i <= ki; ++i) {
                double temp1 = i;
                w(i + 1, 1) = 1.0 / temp1;
            }

            // Compute g[*]
            double term = 0.0;
            for (int j = 2; j <= ki; ++j) {
                double psijm1 = psi_(j, 1);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 1; i <= ki + 1 - j; ++i) {
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                }
                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }

            // Interpolate for the solution yout and for
            // the derivative of the solution ypout
            for (int j = 1; j <= ki; ++j) {
                int i = ki + 1 - j;
                for (int idx = 0; idx < n_eqn; ++idx) {
                    yout(idx, 1) += g(i + 1, 1) * phi(idx, i + 1);
                    ypout(idx, 1) += rho(i + 1, 1) * phi(idx, i + 1);
                }
            }

            yout = y + hi * yout;
            y = yout;
            State_ = DE_STATE::DE_DONE; // % Set return code
            t = tout; // Set independent variable
            told = t; // Store independent variable
            OldPermit = PermitTOUT;
            return y; // Normal exit
        }























        //% If cannot go past output point and sufficiently close,
        //% extrapolate and return
        if (~PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
            double h = tout - x;
            yp = func(x, yy);          //% Compute
            double derivative = yp(x, 1);
            y = yy + h * yp;                //% Extrapolate vector from x to tout
            State_ = DE_STATE::DE_DONE; //% Set return code
            t = tout;             //% Set independent variable
            told = t;                //% Store independent variable
            OldPermit = PermitTOUT;
            return y;                       //% Normal exit
        }

//       % Test for too much work
//       %   if (nostep >= maxnum)
//       %       State_ = DE_STATE.DE_NUMSTEPS; % Too many steps
//       %       if (stiff)
//       %           State_ = DE_STATE.DE_STIFF;% Stiffness suspected
//       %       end
//       %       y         = yy;                % Copy last step
//                                                                                                             %       t         = x;
//       %       told      = t;
//       %       OldPermit = true;
//       %       return;                        % Weak failure exit
//       %   end

        //% Limit step size, set weight vector and take a step
        double h = sign_(fmin(fabs(h), fabs(tend - x)), h);
        for (int l = 1; l <= n_eqn; l++) {
            wt(l, 1) = releps * fabs(yy(l, 1)) + abseps;
        }

        //%   Step
        //%
        //% Begin block 0
        //%
        //% Check if step size or error tolerance is too small for machine
        //% precision.  If first step, initialize phi array and estimate a
        //% starting step size. If step size is too small, determine an
        //% acceptable one.


        if (fabs(h) < fouru * fabs(x)) {
            h = sign_(fouru * fabs(x), h);
            bool crash = true;
            return y;           //% Exit
        }

        double p5eps = 0.5 * epsilon;
        bool crash = false;
        g(2, 1) = 1.0;
        g(3, 1) = 0.5;
        sig(2, 1) = 1.0;

        int ifail = 0;

        //% If error tolerance is too small, increase it to an
        //% acceptable value.

        double round = 0.0;
        for (int l = 1; l <= n_eqn; l++) {
            round = round + (y(l, 1) * y(l, 1)) / (wt(l, 1) * wt(l, 1));
        }
        round = twou * sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            return y;
        }

        bool nornd;
        int k;
        bool phase1;
        double absh;
        double hnew;

        if (start) {
            //% Initialize.Compute appropriate step size for first step.
            yp = func(x, y);
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum = sum + (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));
            }
            sum = sqrt(sum);
            absh = fabs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * sqrt(epsilon / sum);
            }
            h = sign_(fmax(absh, fouru * fabs(x)), h);
            double hold = 0.0;
            hnew = 0.0;
            k = 1;
            int kold = 0;
            start = false;
            phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; l++) {
                    phi(l, 16) = 0.0;
                }
            }
        }


        //% End block 0


        //% Repeat blocks 1, 2 (and 3) until step is successful

        int kp1;
        int kp2;
        int km1;
        int km2;
        int knew;
        double erk;
        double erkm1;
        double erkm2;
        int ns;
        double temp2;

        while (true) {


            //% Begin block 1
            //% Compute coefficients of formulas for this step.Avoid computing
            //% those quantities not changed when step size is not changed.
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;
            //% ns is the number of steps taken with size h, including the
            //% current one. When k<ns, no coefficients change.

            int hold;
            if (h != hold) {
                ns = 0;
            }
            int kold;
            if (ns <= kold) {
                ns = ns + 1;
            }
            int nsp1 = ns + 1;


            if (k >= ns) {
                //% Compute those components of alpha[*],beta[*],psi[*],sig[*]
                //% which are changed
                beta(ns + 1, 1) = 1.0;
                int realns = ns;
                alpha(ns + 1, 1) = 1.0 / realns;
                double temp1 = h * realns;
                sig(nsp1 + 1, 1) = 1.0;
                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; i++) {
                        int im1 = i - 1;
                        temp2 = psi_(im1 + 1, 1);
                        psi_(im1 + 1, 1) = temp1;
                        beta(i + 1, 1) = beta(im1 + 1, 1) * psi_(im1 + 1, 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i + 1, 1) = h / temp1;
                        int reali = i;
                        sig(i + 2, 1) = reali * alpha(i + 1, 1) * sig(i + 1, 1);
                    }
                }
                psi_(k + 1, 1) = temp1;

                //% Compute coefficients g[*]; initialize v[*] and set w[*].

                if (ns > 1) {
                    //% If order was raised, update diagonal part of v[*]
                    if (k > kold) {
                        int temp4 = k * kp1;
                        v(k + 1, 1) = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 1; j <= nsm2; j++) {
                            int i = k - j;
                            v(i + 1, 1) = v(i + 1, 1) - alpha(j + 2, 1) * v(i + 2, 1);
                        }
                    }

                    //% Update V[*] and set W[*]
                    int limit1 = kp1 - ns;
                    double temp5 = alpha(ns + 1, 1);
                    for (int iq = 1; iq <= limit1; iq++) {
                        v(iq + 1, 1) = v(iq + 1, 1) - temp5 * v(iq + 2, 1);
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }
                    g(nsp1 + 1, 1) = w(2, 1);
                } else {
                    for (int iq = 1; iq <= k; iq++) {
                        int temp3 = iq * (iq + 1);
                        v(iq + 1, 1) = 1.0 / temp3;
                        w(iq + 1, 1) = v(iq + 1, 1);
                    }

                }


                //% Compute the g[*] in the work vector w[*]
                int nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        int limit2 = kp2 - i;
                        double temp6 = alpha(i, 1);
                        for (int iq = 1; iq <= limit2; iq++) {
                            w(iq + 1, 1) = w(iq + 1, 1) - temp6 * w(iq + 2, 1);
                        }
                        g(i + 1, 1) = w(2, 1);
                    }
                }
            }//%if K >= NS


            //% End block 1



            //% Begin block 2
            //% Predict a solution p[*], evaluate derivatives using predicted
            //% solution, estimate local error at order k and errors at orders
            //% k, k-1, k-2 as if constant step size were used.
            //% Change phi to phi star
            if (k >= nsp1) {
                for (int i = nsp1; i <= k; i++) {
                    double temp1 = beta(i + 1, 1);
                    for (int l = 1; l <= n_eqn; l++) {
                        phi(l, i + 1) = temp1 * phi(l, i + 1);
                    }
                }
            }

            //% Predict solution and differences
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, kp2 + 1) = phi(l, kp1 + 1);
                phi(l, kp1 + 1) = 0.0;
                p(l, 1) = 0.0;
            }
            for (int j = 1; j <= k; j++) {
                int i = kp1 - j;
                int ip1 = i + 1;
                temp2 = g(i + 1, 1);
                for (int l = 1; l <= n_eqn; l++) {
                    p(l, 1) = p(l, 1) + temp2 * phi(l, i + 1);
                    phi(l, i + 1) = phi(l, i + 1) + phi(l, ip1 + 1);
                }
            }
            if (nornd) {
                p = y + h * p;
            } else {
                for (int l = 1; l <= n_eqn; l++) {
                    double tau = h * p(l, 1) - phi(l, 16);
                    p(l, 1) = y(l, 1) + tau;
                    phi(l, 17) = (p(l, 1) - y(l, 1)) - tau;
                }
            }
            double xold = x;
            x = x + h;
            absh = fabs(h);
            yp = func(x, p);

            //% Estimate errors at orders k, k-1, k-2
            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;

            for (int l = 1; l <= n_eqn; l++) {
                double temp3 = 1.0 / wt(l, 1);
                double temp4 = yp(l, 1) - phi(l, 1 + 1);
                if (km2 > 0) {
                    erkm2 = erkm2 + ((phi(l, km1 + 1) + temp4) * temp3) * ((phi(l, km1 + 1) + temp4) * temp3);
                }
                if (km2 >= 0) {
                    erkm1 = erkm1 + ((phi(l, k + 1) + temp4) * temp3) * ((phi(l, k + 1) + temp4) * temp3);
                }
                erk = erk + (temp4 * temp3) * (temp4 * temp3);
            }

            if (km2 > 0) {
                erkm2 = absh * sig(km1 + 1, 1) * gstr[km2 + 1] * sqrt(erkm2);
            }
            if (km2 >= 0) {
                erkm1 = absh * sig(k + 1, 1) * gstr[km1 + 1] * sqrt(erkm1);
            }

            double temp5 = absh * sqrt(erk);
            double err = temp5 * (g(k + 1, 1) - g(kp1 + 1, 1));
            erk = temp5 * sig(kp1 + 1, 1) * gstr[k + 1];
            knew = k;

            //% Test if order should be lowered
            if (km2 > 0) {
                if (fmax(erkm1, erkm2) <= erk) {
                    knew = km1;
                }
            }
            if (km2 == 0) {
                if (erkm1 <= 0.5 * erk) {
                    knew = km1;
                }
            }


            //% End block 2
            //% If step is successful continue with block 4, otherwise repeat
            //% blocks 1 and 2 after executing block 3

            bool success = (err <= epsilon);

            if (!success) {

                //% Begin block 3
                //% The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
                //% 3rd consecutive failure, set order to 1. If step fails more
                //% than 3 times, consider an optimal step size. Double error
                //% tolerance and return if estimated step size is too small
                //% for machine precision.
                //% Restore x, phi[*,*] and psi[*]
                phase1 = false;
                x = xold;
                for (int i = 1; i <= k; i++) {
                    double temp1 = 1.0 / beta(i + 1, 1);
                    int ip1 = i + 1;
                    for (int l = 1; l <= n_eqn; l++) {
                        phi(l, i + 1) = temp1 * (phi(l, i + 1) - phi(l, ip1 + 1));
                    }
                }

                if (k >= 2) {
                    for (int i = 2; i <= k; i++) {
                        psi_(i, 1) = psi_(i + 1, 1) - h;
                    }
                }

                //% On third failure, set order to one.
                //% Thereafter, use optimal step size
                ifail = ifail + 1;
                temp2 = 0.5;
                if (ifail > 3) {
                    if (p5eps < 0.25 * erk) {}
                    temp2 = sqrt(p5eps / erk);
                }

                if (ifail >= 3) {
                    knew = 1;
                }
                h = temp2 * h;
                k = knew;
                if (fabs(h) < fouru * fabs(x)) {
                    crash = true;
                    h = sign_(fouru * fabs(x), h);
                    epsilon = epsilon * 2.0;
                    return y;                 //% Exit
                }


                //% End block 3, return to start of block 1

            }//% end if (success)

            if (success) {
                break;
            }

        }


        //% Begin block 4
        //% The step is successful. Correct the predicted solution, evaluate
        //% the derivatives using the corrected solution and update the
        //% differences. Determine best order and step size for next step.

        double kold = k;
        double hold = h;

        //% Correct and evaluate
        double temp1 = h * g(kp1 + 1, 1);
        if (nornd) {
            for (int l = 1; l <= n_eqn; l++) {
                y(l, 1) = p(l, 1) + temp1 * (yp(l, 1) - phi(l, 2));
            }
        } else {
            for (int l = 1; l <= n_eqn; l++) {
                rho(l, 1) = temp1 * (yp(l, 1) - phi(l, 2)) - phi(l, 17);
                y(l, 1) = p(l, 1) + rho(l, 1);
                phi(l, 16) = (y(l, 1) - p(l, 1)) - rho(l, 1);
            }
        }
        yp = func(x, y);

        //% Update differences for next step
        for (int l = 1; l <= n_eqn; l++) {
            phi(l, kp1 + 1) = yp(l, 1) - phi(l, 2);
            phi(l, kp2 + 1) = phi(l, kp1 + 1) - phi(l, kp2 + 1);
        }
        for (int i = 1; i <= k; i++) {
            for (int l = 1; l <= n_eqn; l++) {
                phi(l, i + 1) = phi(l, i + 1) + phi(l, kp1 + 1);
            }
        }

        //% Estimate error at order k+1 unless
        //% - in first phase when always raise order,
        //% - already decided to lower order,
        //% - step size not constant so estimate unreliable
        double erkp1 = 0.0;
        if ((knew == km1) || (k == 12)) {
            phase1 = false;
        }

        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else {
            if (knew == km1) {
                //% lower order
                k = km1;
                erk = erkm1;
            } else {
                if (kp1 <= ns) {
                    for (int l = 1; l <= n_eqn; l++) {
                        erkp1 = erkp1 + (phi(l, kp2 + 1) / wt(l, 1)) * (phi(l, kp2 + 1) / wt(l, 1));
                    }
                    erkp1 = absh * gstr[kp1 + 1] * sqrt(erkp1);
                    //% Using estimated error at order k+1, determine
                    //% appropriate order for next step
                    if (k > 1) {
                        if (erkm1 <= fmin(erk, erkp1)) {
                            //% lower order
                            k = km1;
                            erk = erkm1;
                        } else {
                            if ((erkp1 < erk) && (k != 12)) {
                                //% raise order
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    } else if (erkp1 < 0.5 * erk) {
                        //% raise order
                        //% Here erkp1 < erk < max(erkm1,ermk2) else
                        //% order would have been lowered in block 2.
                        //% Thus order is to be raised
                        k = kp1;
                        erk = erkp1;
                    }
                }//% end if kp1 <= ns
            }//% end if knew != km1
        } //% end if !phase1

        //% With new order determine appropriate step size for next step
        if (phase1 || (p5eps >= erk * two[k + 2])) {
            hnew = 2.0 * h;
        } else {
            if (p5eps < erk) {
                temp2 = k + 1;
                double r = p5eps / pow(erk, 1.0 / temp2);
                hnew = absh * fmax(0.5, fmin(0.9, r));
                hnew = sign_(fmax(hnew, fouru * fabs(x)), h);
            } else {
                hnew = h;
            }
        }
        h = hnew;


        //% End block 4
        //% Test for too small tolerances
        if (crash) {
            State_ = DE_STATE::DE_BADACC;
            relerr = epsilon * releps;       //% Modify relative andabsolute
            abserr = epsilon * abseps;       //% accuracy requirements
            y = yy;                   //% Copy last step
            t = x;
            told = t;
            OldPermit = true;
            return y;                       //% Weak failure exit
        }

        nostep = nostep + 1;  //% Count total number of steps
        //% Count number of consecutive steps taken with the order of
        kle4 = kle4 + 1;
        if (kold > 4) {
            kle4 = 0;
        }
        if (kle4 >= 50) {
            stiff = true;
        }


    }


        //% End step loop
        //%   if ( State_==DE_STATE.DE_INVPARAM )
        //%       error ('invalid parameters in DEInteg');
        //%       exit;
        //%   end
        //%   if ( State_==DE_STATE.DE_BADACC )
        //%       warning ('on','Accuracy requirement not achieved in DEInteg');
        //%   end
        //%   if ( State_==DE_STATE.DE_STIFF )
        //%       warning ('on','Stiff problem suspected in DEInteg');
        //%   end
        //%   if ( State_ >= DE_STATE.DE_DONE )
        //%       break;
        //%   end
        //%
        //% end



}