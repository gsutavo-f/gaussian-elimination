use std::time::Instant;

fn main() {
    let n: usize = 3;
    let a = vec![
        vec![2.0, 1.0, -1.0],
        vec![1.0, 2.0, 1.0],
        vec![1.0, 1.0, 1.0],
    ];
    let b = vec![-3.0, 3.0, 2.0];
    let x = vec![0.0; 3];
    _gauss_solver_execution_time(n, a, b, x);

    let n: usize = 3;
    let a = vec![
        vec![1.0, 1.0, 1.0],
        vec![-2.0, 1.0, 1.0],
        vec![1.0, 3.0, 1.0],
    ];
    let b = vec![2.0, 5.0, 4.0];
    let x = vec![0.0; 3];
    _gauss_solver_execution_time(n, a, b, x);

    let n: usize = 3;
    let a = vec![
        vec![3.0, 2.0, -1.0],
        vec![2.0, -2.0, 4.0],
        vec![-1.0, 0.5, -1.0],
    ];
    let b = vec![1.0, -2.0, 0.0];
    let x = vec![0.0; 3];
    _gauss_solver_execution_time(n, a, b, x);

    let n: usize = 2;
    let a = vec![vec![2.0, 3.0], vec![-3.0, 3.0]];
    let b = vec![6.0, 15.0];
    let x = vec![0.0; 2];
    _gauss_solver_execution_time(n, a, b, x);

    let n: usize = 4;
    let a = vec![
        vec![4.0, 1.0, 2.0, -3.0],
        vec![-3.0, 3.0, -1.0, 4.0],
        vec![-1.0, 2.0, 5.0, 1.0],
        vec![5.0, 4.0, 3.0, -1.0],
    ];
    let b = vec![-16.0, 20.0, -4.0, -10.0];
    let x = vec![0.0; 4];
    _gauss_solver_execution_time(n, a, b, x);

    let n: usize = 5;
    let a = vec![
        vec![4.0, 1.0, 2.0, -3.0, 5.0],
        vec![-3.0, 3.0, -1.0, 4.0, -2.0],
        vec![-1.0, 2.0, 5.0, 1.0, 3.0],
        vec![5.0, 4.0, 3.0, -1.0, 2.0],
        vec![1.0, -2.0, 3.0, -4.0, 5.0],
    ];
    let b = vec![-16.0, 20.0, -4.0, -10.0, 3.0];
    let x = vec![0.0; 5];
    _gauss_solver_execution_time(n, a, b, x);
}

fn _gauss_solver_execution_time(n: usize, a: Vec<Vec<f64>>, b: Vec<f64>, x: Vec<f64>) {
    let start = Instant::now();
    gauss_solver(n, a, b, x);
    let duration = start.elapsed();

    println!("Time elapsed in gauss_solver() is: {:?}\n", duration);
}

fn gauss_solver(n: usize, mut a: Vec<Vec<f64>>, mut b: Vec<f64>, mut x: Vec<f64>) {
    //ETAPA DE ESCALONAMENTO
    for k in 0..n - 1 {
        let mut max: f64 = a[k][k].abs();
        let mut max_index: usize = k;

        for i in k + 1..n {
            if max < a[i][k].abs() {
                max = a[i][k].abs();
                max_index = i;
            }
        }

        if max_index != k {
            for j in 0..n {
                let temp: f64 = a[k][j];
                a[k][j] = a[max_index][j];
                a[max_index][j] = temp;
            }
            let temp: f64 = b[k];
            b[k] = b[max_index];
            b[max_index] = temp;
        }

        if a[k][k] == 0.0 {
            println!("A matriz dos coeficientes é singular\n");
            return;
        } else {
            for m in k + 1..n {
                let f: f64 = -a[m][k] / a[k][k];
                a[m][k] = 0.0;
                b[m] = b[m] + f * b[k];
                for l in k + 1..n {
                    a[m][l] = a[m][l] + f * a[k][l];
                }
            }
        }
    }

    //ETAPA DE RESOLUÇÃO DO SISTEMA
    for i in (0..=n - 1).rev() {
        x[i] = b[i];
        for j in i + 1..n {
            x[i] = x[i] - x[j] * a[i][j];
        }
        x[i] = x[i] / a[i][i];
    }

    //IMPRIME RESULTADO
    for (i, value) in x.iter().enumerate() {
        println!("x{} = {}", i, value);
    }
}
