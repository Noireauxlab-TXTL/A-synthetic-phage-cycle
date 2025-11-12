for (i = 1; i <= 6; i++) {
    camera_position = "s" + i + "_"; // String concatenation

    // Open and process gray channel
    run("Image Sequence...", "open=[/Volumes/LaCie/absolute VEZ INFECTION/08-24-24] file=[BF Camera_" + camera_position + "] sort use");
    saveAs("Tiff", "/Volumes/LaCie/absolute VEZ INFECTION/08-24-24/films/ziane A1D3/s" + camera_position + "bf.tif");

    run("Image Sequence...", "open=[/Volumes/LaCie/absolute VEZ INFECTION/08-24-24] file=[GFP Camera_" + camera_position + "] sort use");
    saveAs("Tiff", "/Volumes/LaCie/absolute VEZ INFECTION/08-24-24/films/ziane A1D3/s" + camera_position + "gfp.tif");

    
    //run("Image Sequence...", "open=[/Volumes/LaCie/absolute VEZ INFECTION/12-12-24] file=[BF Camera_" + camera_position + "] sort use");
    //run("Image Sequence...", "open=/Volumes/LaCie/absolute VEZ INFECTION/11-29-24/edta ves infect lysis/", "virtual filter=[BF Camera_" + camera_position + "] start=2 count=26");

    // Open and process red channel
    run("Image Sequence...", "open=[/Volumes/LaCie/absolute VEZ INFECTION/08-24-24] file=[TexRed Camera_" + camera_position + "] sort use");
    saveAs("Tiff", "/Volumes/LaCie/absolute VEZ INFECTION/08-24-24/films/ziane A1D3/s" + camera_position + "red.tif");
    //run("Image Sequence...", "open=/Volumes/LaCie/absolute VEZ INFECTION/11-29-24/edta ves infect lysis/", "virtual filter=[TexRed Camera_" + camera_position + "] start=2 count=26");
}