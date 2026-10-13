import java.util.Scanner;

class Complex{
    private final double re;
    private final double im;

    public Complex(double re, double im) {
        this.re = re;
        this.im = im;
    }

    public double getRe() { return re; }
    public double getIm() { return im; }

    public Complex add(Complex o) {
        return new Complex(re + o.re, im + o.im);
    }

    public Complex subtract(Complex o) {
        return new Complex(re - o.re, im - o.im);
    }

    public Complex multiply(Complex o) {
        return new Complex(re * o.re - im * o.im,
                re * o.im + im * o.re);
    }

    public Complex divide(Complex o) {
        double d = o.re * o.re + o.im * o.im;
        if (d == 0) throw new ArithmeticException("Деление на ноль");
        return new Complex((re * o.re + im * o.im) / d,
                (im * o.re - re * o.im) / d);
    }

    public boolean isZero() {
        return Math.abs(re) < 1e-12 && Math.abs(im) < 1e-12;
    }

    @Override
    public String toString() {
        if (im == 0) return String.valueOf(re);
        if (re == 0) return im + "i";
        return re + (im > 0 ? "+" : "-") + Math.abs(im) + "i";
    }
}

class ComplexMatrix{
    private final int x;
    private final int y;
    private Complex[][] data;


    public ComplexMatrix(int x,int y){
        this.x=x;
        this.y=y;
        this.data= new Complex[x][y];
        for (int i=0;i<x;++i){
            for (int j=0;j<y;++j){
                data[i][j]=new Complex(0,0);
            }
        }
    }
    public Complex getMat(int i,int j){return data[i][j];}
    public void setMat(int i,int j,Complex value){data[i][j]=value;}


    public ComplexMatrix SumMatrix (ComplexMatrix mat1){
        ComplexMatrix sum=new ComplexMatrix(this.x,this.y);
        if (mat1.x!=this.x || mat1.y!=this.y ) throw new IllegalArgumentException("Матрицы разного размера");
        for (int i=0; i<this.x;++i){
            for (int j=0; j<this.y;++j){
                sum.data[i][j]= mat1.data[i][j].add(this.data[i][j]) ;
            }
        }
        return sum;
    }


    public ComplexMatrix CompMartix(ComplexMatrix mat1){
        ComplexMatrix proizv=new ComplexMatrix(this.x,mat1.y);
        if (mat1.x != this.y) throw new IllegalArgumentException("Матрицы разного размера");

        for (int i=0; i<this.x;++i){
            for (int j=0; j<mat1.y;++j){
                Complex s=new Complex(0,0);
                for (int l=0;l<this.y;++l){
                   s=s.add(this.data[i][l].multiply(mat1.data[l][j]));
                }
                proizv.data[i][j]=s;
            }
        }
        return proizv;
    }
    Complex Determinant () {
        if (x != y) throw new IllegalArgumentException("Определитель есть только у квадратной матрицы");
        int n=x;

        Complex[][] a= new Complex[n][n];
        for (int i = 0; i<n;++i){
            for(int j=0;j<n;++j){
                a[i][j]=data[i][j];
            }
        }
        Complex det=new Complex(1,0  );
        for (int k = 0; k < n; ++k) {
            int p=k;
            while (p<n && a[p][k].isZero()) ++p;
            if (p==n)return new Complex(0,0);
            if (p!=k) {
                Complex[]tmp=a[p];
                a[p]=a[k];
                a[k]=tmp;
                det=det.multiply(new Complex(-1,0));

            }

            for (int i = k + 1; i < n; ++i) {
                Complex coef = a[i][k].divide(a[k][k]);
                for (int j = k; j < n; ++j) {
                    a[i][j] = a[i][j].subtract(coef.multiply(a[k][j]));
                }
            }
            det = det.multiply(a[k][k]);
        }
        return det;
        }


    public ComplexMatrix inverse() {
        if (x != y)
            throw new IllegalArgumentException("Обратная матрица есть только у квадратной");

        int n = x;


        Complex[][] a = new Complex[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[i][j] = data[i][j];


        Complex[][] inv = new Complex[n][n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                inv[i][j] = (i == j) ? new Complex(1, 0) : new Complex(0, 0);

        for (int k = 0; k < n; ++k) {

            int p = k;
            while (p < n && a[p][k].isZero()) ++p;
            if (p == n)
                throw new ArithmeticException("Определитель равен 0, обратной матрицы нет");

            if (p != k) {
                Complex[] tmp = a[p]; a[p] = a[k]; a[k] = tmp;
                tmp = inv[p]; inv[p] = inv[k]; inv[k] = tmp;
            }


            Complex lead = a[k][k];
            for (int j = 0; j < n; ++j) {
                a[k][j] = a[k][j].divide(lead);
                inv[k][j] = inv[k][j].divide(lead);
            }


            for (int i = 0; i < n; ++i) {
                if (i == k) continue;
                Complex coef = a[i][k];
                for (int j = 0; j < n; ++j) {
                    a[i][j] = a[i][j].subtract(coef.multiply(a[k][j]));
                    inv[i][j] = inv[i][j].subtract(coef.multiply(inv[k][j]));
                }
            }
        }

        ComplexMatrix result = new ComplexMatrix(n, n);
        result.data = inv;
        return result;
    }

    public ComplexMatrix DelMatr(ComplexMatrix B) {
        return this.CompMartix(B.inverse());
    }

    ComplexMatrix Trans() {
        ComplexMatrix result = new ComplexMatrix(this.y, this.x);
        for (int i = 0; i < this.x; ++i)
            for (int j = 0; j < this.y; ++j)
                result.data[j][i] = this.data[i][j];
        return result;
    }
    @Override

    public String toString(){
        StringBuilder st =new StringBuilder();
        for(int i =0;i<x;++i){
            for(int j =0;j<y;++j){
                st.append(data[i][j]);
                st.append("\t");
            }
            st.append("\n");
        }
        return st.toString();
    }


}
public class Main{

    static ComplexMatrix readMatrix(Scanner sc){
        System.out.print("Число строк:");
        int rows = sc.nextInt();
        System.out.print("Число столбцов:");
        int cols = sc.nextInt();
        if (rows <= 0 || cols <= 0)
            throw new IllegalArgumentException("Размеры должны быть больше нуля");

        ComplexMatrix m = new ComplexMatrix(rows,cols);
        System.out.println("Вводите элементы построчно, каждый как два числа: re im");
        for (int i = 0; i<rows;++i) {
            for (int j = 0; j<cols;++j) {
                System.out.print("[" + i + "][" + j + "] = ");
                double re = sc.nextDouble();
                double im = sc.nextDouble();
                m.setMat(i,j,new Complex(re,im));
            }
        }
        return m;
    }

    public static void main(String[] args){
        Scanner sc=new Scanner(System.in);
        System.out.println("Введите матрицу A");
        ComplexMatrix A=readMatrix(sc);

        while (true){
            System.out.println("\n1 - сложить A + B");
            System.out.println("2 - умножить A * B");
            System.out.println("3 - разделить A / B");
            System.out.println("4 - транспонировать A");
            System.out.println("5 - определитель A");
            System.out.println("0 - выход");
            System.out.print("Выбор: ");
            int choice = sc.nextInt();
            try{
                switch (choice){
                    case 1->{
                        System.out.println("Введите матрицу B");
                        ComplexMatrix B=readMatrix(sc);
                        System.out.println("Результат:\n"+A.SumMatrix(B));
                    }
                    case 2->{
                        System.out.println("Введите матрицу B");
                        ComplexMatrix B=readMatrix(sc);
                        System.out.println("Результат:\n"+A.CompMartix(B));
                    }
                    case 3->{
                        System.out.println("Введите матрицу B");
                        ComplexMatrix B=readMatrix(sc);
                        System.out.println("Результат:\n"+A.DelMatr(B));
                    }
                    case 4->{
                        System.out.println(A.Trans());
                    }
                    case 5->{
                        System.out.println("det = " + A.Determinant());
                    }
                    case 0-> {return;}
                    default -> System.out.println("Нет такого пункта");
                }
            }catch (IllegalArgumentException|ArithmeticException e){
                System.out.println("Ошибка: "+ e.getMessage());
            }
        }
    }
}
