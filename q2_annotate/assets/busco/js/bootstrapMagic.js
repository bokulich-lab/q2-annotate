function removeBS3refs() {
    // remove Bootstrap 3 CSS/JS reference
    let head = document.getElementsByTagName("head")[0]
    let links = head.getElementsByTagName("link")
    for (let i = 0; i < links.length; i++) {
        if (links[i].href.includes("q2templateassets/css/bootstrap")) {
            links[i].remove()
        }
    }
    let scripts = head.getElementsByTagName("script")
    for (let i = 0; i < scripts.length; i++) {
        if (scripts[i].src.includes("q2templateassets/js/bootstrap")) {
            scripts[i].remove()
        }
    }
}

function adjustTagsToBS3() {
    // adjust tags to BS3
    let tabs = document.getElementsByClassName("nav nav-tabs")[0].children
    for (let i = 0; i < tabs.length; i++) {
        let isActive = tabs[i].className.includes("active")
        tabs[i].className = "nav-item"
        let link = tabs[i].getElementsByTagName("a")[0]
        if (isActive) {
            link.classList.add("active")
        }
        link.classList.add("nav-link")

    }
}
