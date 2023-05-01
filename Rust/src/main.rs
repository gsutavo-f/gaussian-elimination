fn main() {
    let n :i32 = 3;
    let a :[[i32;3];3] = [[2, 1, -1], [1, 2, 1], [1, 1, 1]];
    let b :[i32;3] = [-3, 3, 2];
    let x :[i32;3] = [0;3];

    gauss_solver(n, a, b, x);
}

fn gauss_solver(n :i32, a :[[i32;3];3], b :[i32;3], x :[i32;3]) {
    println!("n={n}");

    for (i, row) in a.iter().enumerate() {
        for (j, col) in row.iter().enumerate() {
            println!("[{}][{}]={}", i, j, col);
        }
    }

    for (i, row) in b.iter().enumerate() {
        println!("[{}]={}", i, row);
    }

    for (i, row) in x.iter().enumerate() {
        println!("[{}]={}", i, row);
    }
}
