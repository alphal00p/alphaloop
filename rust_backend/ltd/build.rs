fn main() {
    println!("cargo:rustc-link-search=../libraries/fjcore");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");
}