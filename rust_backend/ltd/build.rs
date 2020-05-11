fn main() {
    println!("cargo:rustc-link-search=../libraries/fjcore");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");

    #[cfg(feature = "mg_numerator")]
    {
        let path = std::env::var("MG_NUMERATOR_PATH")
            .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

        println!("cargo:rustc-link-search={}/lib", path);
        println!("cargo:rustc-link-lib=dhelas");
        println!("cargo:rustc-link-lib=model");
        println!("cargo:rustc-link-lib=MGnumerators");
        println!("cargo:rustc-link-lib=gfortran");
    }
}
