//Vallejos Bardales Gian - June 2014
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <conio.h>
#define REP(i, a) for( int i = 0; i < a; i++ )
#define RFOR(i,x,y) for(int i = x; i>= y; i--)
#define FOR(i,x,y) for (int i = x; i <= y; i++)
#define ITFOR(it,A) for(typeof A.begin() it = A.begin(); it!=A.end(); it++)
#define all(v) (v).begin(), (v).end()
#define rall(v) (v).rbegin(), (v).rend()
#define VE vector <int>
#define mset(A,x) memset(A, x, sizeof A)
#define PB push_back
#define ones(x) __builtin_popcount(x)
#define cua(x) (x)*(x)
#define debug(x) cout <<#x << " = " << x << endl
#define adebug(x,n) cout <<#x<<endl; REP(i,n)cout<<x[i]<<char(i+1==n?10:32)
#define mdebug(x,m,n) cout <<#x<<endl; REP(i,m)REP(j,n)cout<<x[i][j]<<char(j+1==n?10:32)
#define to_string(A) static_cast<ostringstream*>(&(ostringstream() << A))->str()
using namespace std;
#define ll long long
#define N 20005
#define B 100
//VARIABLES GLOBALES
vector <ll> factorBase, numLeft;
const int _size = 10000;
ll matriz[_size+2][_size+2];
ll inversa[_size+2][_size+2];
ll n_row, m_col, posAnsSize;
vector <string> ans;
vector <ll> posAns[_size+2];
ll n, x, y;
//CRIBA ERATOSTENES
template <class T> class CribaEratostenes{
public:
    vector <bool> primos;
    void aplicarCribaEratostenes(T n){
        primos.assign(n+10, false);
        FOR( i, 2, sqrt(n) ){
            if( !primos[i] )
                FOR( j, i, n/i )  primos[i*j] = true;
        }
    }
};
template <class T> inline T gcd( T a, T b ){ return (!b) ? a : gcd(b, a%b); }
//JACOBI
inline ll getExp( ll a ){
    int cnt = 0;
    while( a%2 == 0 ){ a /= 2; cnt++; }
    return cnt;
}
ll jacobi( ll a, ll n ){
    if( a == 0 || a == 1)
        return a;
    ll a1 = a / pow(2, getExp(a));
    ll s, n1;
    if( getExp(a) % 2 == 0 ){
        s = 1;
    }else{
        if( (n-1)%8 == 0 || (n-7)%8 == 0 ){ s = 1; }
		else if( (n-3)%8 == 0 || (n-5)%8 == 0 ){  s = -1; }
    }
    if( (n-3)%4 == 0 && (a1-3)%4 == 0 ){  s = s * -1;  }
    n1 = n%a1;
    return ( a1 == 1 ) ? s : (s*jacobi(n1,a1));
}
//MOSTRAR MATRIZ MOD 2
void mostrarMatrizMod2(){
    cout << "FB: ";
    REP(i, factorBase.size() )
        cout << factorBase[i] << " ";
    cout << endl;
    REP(i, numLeft.size() ){
        cout << numLeft[i] << ": ";
        REP(j, factorBase.size()){
            cout << matriz[i][j] << " ";
        }
        cout << endl;
    }
}
//FACTOR BASE
void generarFactorBase(){
    factorBase.PB(-1); //Factor Base inicial -1
    //Criba Eratostenes
    CribaEratostenes <ll> ce;
    ce.aplicarCribaEratostenes(B);
  //Recorrer B valores
    FOR(i, 2, B){
        if( !ce.primos[i] && jacobi(n, i) == 1 ) factorBase.PB(i); //Si es primo y Jacobi se cuenta como Factor base
        //if( factorBase.size() == 8 ) break;
    }
}

bool esFunction(ll n, ll x, ll Fx, int row){
    int exp[factorBase.size()+2]; //Saca exponentes
    FOR(i, 1, factorBase.size() - 1 ){
        exp[i] = 0;
        while( x%factorBase[i] == 0 ){ //Va viendo si esta compuesto por los factores base
            x /= factorBase[i];
            exp[i]++; //Almacena exponentes válidos
        }
    }
    if( x == 1 ){
        matriz[row][0] = ( Fx*Fx < n ) ? 1 : 0; //Si el factorBase -1 existe
        FOR(i, 1, factorBase.size() - 1) //Crear la matriz de acuerdo a los exponentes
            matriz[row][i] = (exp[i]%2 != 0) ? 1 : 0; //Asignar 1 o 0 a la matriz
        return true;
    }else{ return false; }
}

void generarMatriz(){
    ll x, Fx;
    int raiz = sqrt(n), row = 0;
    numLeft.clear(); //Números de la izquierda clear
    FOR(i, 1, B){ //Ve números del 1-B en +-i sqrt(n)
        //if(numLeft.size() == factorBase.size() ) break;
      //raiz - i
        x = raiz - i;
        Fx = (x * x) - n; //Calcular F(n) = n^2 - n
        if(Fx < 0) Fx = Fx * -1; //Si es negativo, pásarle positivo
        if( esFunction(n, Fx, x, row) ){ numLeft.PB(x); row++; }
      //raiz + i
        x = raiz + i;
        Fx = (x * x) - n;
        if(Fx < 0) Fx = Fx * -1;
        if( esFunction(n, Fx, x, row) ){ numLeft.PB(x); row++; }
    }
}

void mostrarMatriz(int X[100+2][100+2], int n, int m){
	REP(i, n){
		REP(j, m){ cout << X[i][j] << " ";	} cout << endl;
	}
}

void generarInversa(int n){
	memset(inversa, 0, sizeof(inversa));
	REP(i, n){ inversa[i][i] = 1; }
}

void sumarMatrizInversa(int f1, int f2){
	for( int col = 0; col < m_col; col++ ){
		matriz[f2][col] = (matriz[f2][col] + matriz[f1][col])%2;
		inversa[f2][col] = (inversa[f2][col] + inversa[f1][col])%2;
	}
}

void gaussModuloDos(int n, int m){
	int posCol;
	for( int row = 0; row < n; row++ ){
		posCol = -1;
		for( int col = 0; col < m; col++ ){
			if( matriz[row][col] == 1 ){ posCol = col; break; }
		}
		if( posCol != -1 ){
			for( int rAux = row + 1; rAux < n; rAux++ ){
				if( matriz[rAux][posCol] == 1 ){ sumarMatrizInversa(row, rAux);	}
			}
		}
	}
}

bool esCero(int row, int m){
	for( int col = 0; col < m; col++ ){
		if( matriz[row][col] == 1 )	return false;
	}
	return true;
}

void guardarSolucion(int row, int m){
	string res = "";
	for( int col = 0; col < n_row; col++ ){
		res += char(inversa[row][col] + 48);
	}
	ans.push_back(res);
}

void evaluarSolucion(int n, int m){
	for( int row = 0; row < n; row++ ){
		if( esCero(row, m) ){ guardarSolucion(row, m); }
	}
}

void mostrarAnswer(){
    // sort(ans.begin(), ans.end());
	for(int i = 0; i < (int)ans.size(); i++ ){	cout << ans[i] << endl;	}
}

void mostrarPosAns(){
    REP(i, ans.size()){
        REP(j, posAns[i].size() ){ cout << posAns[i][j] << " "; } cout << endl;
    }
}

void calcularPosicionesRespuesta(vector <string> ans, int _size){
    REP(i, ans.size()){
        REP(j, ans[i].length() ){ if( ans[i][j] == '1' ) posAns[i].PB(j); }
    }
}

void calcularXY(int row){
    x = 1; y = 1;
    REP(i, posAns[row].size()){
        x *= numLeft[posAns[row][i]];
        y *= ((numLeft[posAns[row][i]]*numLeft[posAns[row][i]])-n);
    }
    y = sqrt(y);
    x = x%n;
    y = y%n;
}

void obtenerResultados(int posAnsSize){
    ll res1, res2;
    bool resultado = false;
    REP(i, posAnsSize){
        calcularXY(i);
        res1 = gcd(x-y, n); res2 = gcd(x+y, n);
        if( res1 * res2 == n && res1 != 1 && res2 != 1){
            cout << "QS of " << n <<": " << res1 << " " << res2 << endl;
            resultado = true;
            break;
        }
    }
    if( !resultado ) printf("There is not solution. Try to change the base.\n");
}

//QUADRATIC SIEVE
void QS(ll a){
    n = a;
  //GENERAR FACTOR BASE
    generarFactorBase();
    generarMatriz();

    //mostrarMatrizMod2();
//GAUSS MODULO 2
    n_row = numLeft.size(); m_col = factorBase.size();
    //cout << "!" << n_row << " " << m_col << endl;
    generarInversa(n_row);
    gaussModuloDos(n_row, m_col);

    evaluarSolucion(n_row, m_col);

   //mostrarAnswer();

//CALCULAR RESPUESTAS
    posAnsSize = ans.size();
    calcularPosicionesRespuesta(ans, posAnsSize);

    obtenerResultados(posAnsSize);
}
void resetPosAns(){
    FOR(i, 0, _size+2){
        posAns[i].clear();
    }
}

void inicializar(){
    memset(matriz, 0, sizeof(matriz));
    memset(inversa, 0, sizeof(inversa));
    ans.clear();
    factorBase.clear();
    numLeft.clear();
    n = x = y = 0;
    n_row = m_col = posAnsSize = 0;
    resetPosAns();
}

int main(){
    freopen("input.txt", "r", stdin);
    //freopen("out.txt", "w", stdout);

    ll number;
    while( ~scanf("%lld", &number)  ){
        inicializar();
        QS(number);
    }

    //fclose(stdin);
    //fclose(stdout);

    return 0;
}
