#pragma once
//Gauss Random

float rand49(void)
{  /*rand_max=7FFF (32767) */
	static int Num = 0;
	double number;
	int    i;
	i = rand();
	number = (double)(i) / ((unsigned)(RAND_MAX));
	Num++;
	if (Num >= RAND_MAX) {
		time_t t;
		t = time(NULL);
		//		srand((unsigned)(t%RAND_MAX));
		Num = 0;
	}
	return (float)number;
}

double Normal(void)
{
	static int iset = 0;
	static double qset;
	double vx, vy, r, temp;
	if (iset == 0)//noise=normal*deviate
	{
		do
		{
			vx = 2.0 * rand49() - 1.0;
			vy = 2.0 * rand49() - 1.0;
			r = vx * vx + vy * vy;
		} while (r >= 1.0 || r == 0);
		temp = sqrt(-2.0 * log(r) / r);
		qset = vy * temp;
		iset = 1;
		//printf("Hello %f\n", vx * temp);
		return (vx * temp);
	}
	else
	{
		iset = 0;
		return qset;
	}
}



double func_Phi(double a)
{
	double result;
	double absolute;
	double tmp;
	absolute = fabs(a);
	if (absolute < pow(10, -16))
		result = -37.42995;
	else if (absolute > 37.43)
		result = -1.110223e-16;
	else
	{
		tmp = (exp(a) - 1) / (exp(a) + 1);
		tmp = fabs(tmp);
		result = log(tmp);
	}
	return result;
}

int func_sign(int x)
{
	int y;
	if (x >= 0)
	{
		y = 1;
	}
	else
	{
		y = -1;
	}
	return y;
}
