// ga.cpp

#include<bits/stdc++.h>

using namespace std;

typedef unsigned long long ull;

const ull inf = ULLONG_MAX;

const int maxngen = 100010;
const int maxn = 100010;

int ngen, n;
int nsur, ncrs, nmut;

//random_device generator;
default_random_engine generator;
uniform_int_distribution<ull> distribution(0,inf);
uniform_int_distribution<int> distbit(0,63);

const float pi = 2 * acosf(0);
const float eps = 1e-7, eps2 = 1e-50;

typedef union {
    ull code;
    struct v {
        float x;
        float y;
    }v;
} chromsome;

vector<chromsome> gen[maxngen];
vector<float> x[maxngen],y[maxngen],f[maxngen], p[maxngen];

float getf(float x, float y) {
    x = fabsf(x);
    y = fabsf(y);
    return (x + y) * (1 + fabs(sinf(x * pi)) + fabs(sinf(y * pi)));
}

bool check(float x, float y) {
    return -60 < x - eps && x + eps < 40 && -30 < y - eps && y + eps < 70;
}

chromsome rand_ch() {
    chromsome t;
    do {
        t.code = distribution(generator);
    } while (!check(t.v.x, t.v.y));
    return t;
}

pair<chromsome, chromsome> crossover(chromsome a, chromsome b) {
    chromsome ta, tb;
    do {
        int i = distbit(generator), j = distbit(generator);
        if ( i > j) swap(i,j);
        ull mask = 1 << i;
        for (int k = i; k < j; ++ k) mask = (mask << 1) + mask;
        mask = (mask & a.code) ^ (mask & b.code);
        ta.code = a.code ^ mask;
        tb.code = b.code ^ mask;
    } while ((!check(ta.v.x, ta.v.y) || !check(tb.v.x, tb.v.y)));
    return {ta, tb};
}

chromsome mutate(chromsome a) {
    chromsome t;
    do {
        int i = distbit(generator);
        ull mask = 1 << i;
        t.code = a.code ^ mask;
    } while (!check(t.v.x, t.v.y));
    return t;
}

int main() {

    ngen = 501;
    n    = 1000;
    nsur = 310;
    ncrs = 680;
    nmut = 10;

    for (int i = 0; i < n; ++ i) {
        gen[0].push_back( rand_ch() );
    }

    for (int g = 0; g < ngen; ++ g) {

        float ans = 1e20, ax, ay;
        float mx = -1e20;
        for (int i = 0; i < n; ++ i) {
            chromsome t = gen[g][i];
            x[g].push_back(t.v.x);
            y[g].push_back(t.v.y);
            f[g].push_back(getf(t.v.x, t.v.y));
            p[g].push_back(log(f[g].back() + eps2));

            if (f[g][i] < ans) {
                ans = f[g][i];
                ax = x[g][i]; ay = y[g][i];
            }
            mx = max(mx, p[g][i]);
        }
        for (int i = 0; i < n; ++ i) {
            p[g][i] = mx + eps - p[g][i];
        }

        if ((g <= 100 && g % 10 == 0) || g % 100 == 0) {
            cout << "[generation " << g << "'s best solution]" << endl;
            cout << "f = " << ans << endl;
            cout << "x = " << ax << endl;
            cout << "y = " << ay << endl << endl;
        }

        discrete_distribution<int> dist(p[g].begin(), p[g].end());
        for (int k = 0; k < nsur; ++ k) {
            int i = dist(generator);
            gen[g + 1].push_back( gen[g][i] );
        }
        for (int k = 0; k < ncrs / 2; ++ k) {
            int i, j;
            i = dist(generator); j = dist(generator);
            auto pr = crossover( gen[g][i], gen[g][j] );
            gen[g + 1].push_back(pr.first);
            gen[g + 1].push_back(pr.second);
        }
        for (int k = 0; k < nmut; ++ k) {
            int i = dist(generator);
            gen[g + 1].push_back( mutate( gen[g][i] ) );
        }

    }

    return 0;
}
